import subprocess
import re
import csv
from collections import defaultdict, Counter

# These three functions
# guess_raw_type, guess_atomic_column_class, guess_column class
# attempt to identify the type of variable in each of the columns
# of a phenotype file

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def guess_raw_type(s):
    # return one of "int","float","str","bool","_empty_"
    if re.match(r'^\s*$', s):
        return "_empty_"
    try:
        f = float(s)
        if "." not in s:
            return "int"
        return "float"
    except ValueError:
        value = s.upper()
        if value =="TRUE" or value == "FALSE":
            return "bool"
        return type(s).__name__

def guess_atomic_column_class(rawtype, obs):
    n_vals = sum(obs.values())
    n_uniq_vals = len(obs)
    if n_vals == n_uniq_vals and rawtype != "float":
        return {"class": "id", "type": rawtype}
    if n_uniq_vals == 1:
        return {"class": "fixed", "type": rawtype, "value": next(iter(obs.keys()))}
    if n_uniq_vals == 2:
        levels = list(obs.keys())
        if rawtype in ["int", "float"]:
            levels.sort(key = lambda x: atof(x))
        return {"class": "binary", "type": rawtype, "levels": levels}
    if rawtype == "str":
        if float(n_uniq_vals)/n_vals > .75:
            return {"class": "descr", "type": "str"}
        else:
            return {"class": "categorical", "type":"str", "levels": list(obs.keys())}
    if rawtype == "float" or rawtype == "int":
        return {"class": "numeric", "type": rawtype}
    return {"class": rawtype, "type": rawtype}

def guess_column_class(colinfo):
    # colinfo is dict (types) of counters (values)
    has_empty = "_empty_" in colinfo and len(colinfo["_empty_"])==1
    col = {k: v for k,v in colinfo.items() if k != "_empty_"}
    n_vals = Counter({k: sum(v.values()) for k, v in col.items()})
    n_uniq_vals = Counter({k: len(v) for k, v in col.items()})
    if len(col)==0:
        #all empty
        return {"class": "fixed", "type": "str"}
    if len(col)==1:
        # all same type
        best_type = n_vals.most_common(1)[0][0]
        ci = guess_atomic_column_class(best_type, col[best_type])
        if has_empty:
            ci["missing"] = next(iter(colinfo["_empty_"].keys()))
        return ci
    if len(col)==2 and "int" in col and "float" in col:
        #promote to float
        best_type = "float"
        vals = col["int"] + col["float"]
        ci =  guess_atomic_column_class(best_type, vals)
        if has_empty:
            ci["missing"] = next(iter(colinfo["_empty_"].keys()))
        return ci
    if len(col)==2 and (n_uniq_vals["str"]==1):
        # likely a single type with a missing indicator
        best_type = [x[0] for x in n_vals.most_common(2) if x[0] != "str"][0]
        ci = guess_atomic_column_class(best_type, colinfo[best_type])
        ci["missing"] = next(iter(colinfo["str"].keys()))
        return ci
    if all((x in ["str","int","float"] for x in col.keys())) and (n_uniq_vals["str"]==1):
        # likely a numeric value with a missing indicator
        best_type = "float"
        vals = col["int"] + col["float"]
        ci = guess_atomic_column_class(best_type, vals)
        ci["missing"] = next(iter(colinfo["str"].keys()))
        return ci
    if len(col)>1 and "str" in col and n_uniq_vals["str"]>2:
        best_type = "str"
        vals = Counter()
        for r in col.values():
            vals += r
        ci = guess_atomic_column_class(best_type, vals)
        if has_empty:
            ci["missing"] = next(iter(colinfo["_empty_"].keys()))
        return ci
    # finally let's just make it a string
    best_type = "str"
    vals = Counter()
    for r in col.values():
        vals += r
    return guess_atomic_column_class(best_type, vals)

def strip_comments(item, token="#"):
    for line in item:
        if not line.strip().startswith(token):
            yield line

def get_comments(item, token="#"):
    for line in item:
        s = line.strip()
        if not s.startswith(token) and len(s)>0:
            raise StopIteration
        yield s

def sniff_file(csvfile):
    chunk = "\n".join([x for _,x in zip(range(50), strip_comments(csvfile))])
    try:
        return csv.Sniffer().sniff(chunk, "\t|, ")
    except:
        return None

def find_header(firstrow, lastcomment, cols):
    colclasses = {k: guess_column_class(v) for k,v in cols.items()}
    col_types = [z["type"] for z in [colclasses[i] for i in range(len(firstrow))]]
    firstrow_types = [guess_raw_type(x) for x in firstrow]

    nonstringcols = sum((x != "str" for x in col_types))
    firstrow_promotions = sum((f=="str" for (f,c) in zip(firstrow_types, col_types) if c!="str"))
    if lastcomment:
        comment_types = [guess_raw_type(x) for x in lastcomment]
        if len(lastcomment) == len(colclasses):
            comment_promotions = sum((f=="str" for (f,c) in zip(comment_types, col_types) if c!="str"))
        else:
            comment_promotions = 0
    else:
        comment_types = []
        comment_promotions = 0

    if comment_promotions > firstrow_promotions and (float)(comment_promotions)/nonstringcols > .9:
        #the last comment has the right number of rows and string in non-string columns
        return (lastcomment, "comment")
    if (float)(firstrow_promotions)/nonstringcols >= .5:
        return (firstrow, "firstrow")
    # no header found, return unique names
    return (["COL{0}".format(i) for i in range(len(colclasses))], "position")

def check_if_ped(cols, obs):
    # This just isn't working quite right
    # better to check column names?
    return False, None
    if len(cols)<6:
        return False, None
    types = [x["type"] for x in cols]
    classes = [x["class"] for x in cols]
    if not (classes[0] == "id" or classes[1]=="id"):
        return False, None
    if not (types[2] == types[1] or classes[2]=="fixed"):
        return False, None
    if not (types[3] == types[1] or classes[3]=="fixed"):
        return False, None
    if not (classes[4]=="binary" or classes[4]=="fixed"):
        return False, None
    return True, None 

def is_possible_sample_id_col(meta, colinfo):
    colclass = meta["class"]
    if meta["class"] == "id" and meta["type"] == "str":
        matches = sum(x.startswith("NWD") for x in colinfo["str"])
        return matches == len(colinfo["str"])
    else:
        return False

def infer_meta(csvfile, dialect=None):
    meta = {"layout": {}, "columns": []}

    start_pos = 0
    csvfile.seek(0)
    first = csvfile.read(4)
    if first.startswith("\xfe\xff") or first.startswith("\xff\xfe"):
        start_pos = 2
    elif first.startswith("\xef\xbb\xbf"):
        start_pos = 3

    # store csv dialect
    if not dialect:
        csvfile.seek(start_pos)
        dialect = sniff_file(csvfile)
    for k in [k for k in dir(dialect) if not k.startswith("_")]:
        meta["layout"]["csv_" + k] = getattr(dialect, k)

    # read comments 
    csvfile.seek(start_pos)
    comments = list(get_comments(csvfile))

    # read and process csv rows 
    csvfile.seek(start_pos)
    cvr = csv.reader(strip_comments(csvfile), dialect)
    firstrow = next(cvr)
    cols = defaultdict(lambda : defaultdict(Counter))
    for row in cvr:
        for idx, val in enumerate(row):
            cols[idx][guess_raw_type(val)][val] += 1
    # find column headers
    if comments:
        lastcomment = next(csv.reader([comments[-1][1:]], dialect))
    else:
        lastcomment = None
    headers, headersource = find_header(firstrow, lastcomment, cols)
    if headersource == "firstrow":
        skip = len(comments)+1
    else:
        skip = len(comments)
    meta["layout"]["skip"] = skip
    # define columns
    meta["columns"] = [None] * len(cols)
    for col, colval in cols.items():
        coldef = guess_column_class(colval)
        coldef["name"] = headers[col]
        meta["columns"][col] = coldef
    #check if ped
    pedlike, ped_columns = check_if_ped(meta["columns"], cols)
    if pedlike:
        #assign standard ped column types
        meta["pedlike"] = 1
        colclasses = ["family_id","sample_id","father_id","mother_id","sex"]
        for actas, col in zip(colclasses, meta["columns"][0:3]):
            if col["class"] != "fixed":
                col["class"] =  actas
    else:
        #if not ped, try to find ID column
        id_col =  next((i for i in range(len(cols)) if  \
            is_possible_sample_id_col(meta["columns"][i], cols[i])), None)
        if id_col is not None:
            meta["columns"][id_col]["class"] = "sample_id"
    return meta

class PhenoReader:
    def __init__(self, path, meta=None):
        self.path = path
        if meta:
            self.meta = meta
        else:
            self.meta = self.infer_meta()

    def infer_meta(self):
        with open(self.path, 'rU') as csvfile:
            return infer_meta(csvfile)

    def get_dialect(self, opts=None):
        class dialect(csv.Dialect):
            pass
        if opts is None and self.meta and "layout" in self.meta:
            opts = {k[4:]:v for (k,v) in self.meta["layout"].items() if k.startswith("csv_")}
            opts = {k: v.encode('utf8') if isinstance(v, str) else v
                for (k,v) in opts.items()} 
            for (k, v) in opts.items():
                setattr(dialect, k, v)
            return dialect
        else:
            return None

    def get_column_missing_values(self):
        cols = self.get_columns();
        return { v["name"]: v.get("missing",None) for v in cols }

    def get_column_indexes(self):
        cols = self.get_columns();
        return { v["name"]: i for (i,v) in enumerate(cols)}

    def get_columns(self):
        if self.meta and self.meta["columns"]:
            return self.meta["columns"]
        else:
            return [];

    def row_extractor(self):
        dialect = self.get_dialect()
        if "layout" in self.meta and "skip" in self.meta["layout"]:
            skip = self.meta['layout']['skip']
        else:
            skip = 0
        with open(self.path, 'rU') as csvfile:
            if not dialect:
                dialect = sniff_file(csvfile)
                csvfile.seek(0)
            if skip>0:
                [csvfile.readline() for i in range(skip)]
            cvr = csv.reader((row for row in csvfile if not row.startswith("#")), dialect)
            for row in cvr:
                yield row

    def get_samples(self):
        sample_col = next(iter([x for x in self.get_columns() if x["class"]=="sample_id"]), None)
        if not sample_col:
            raise "Could not find sample ID column" 
        sample_col_idx = self.get_column_indexes()[sample_col["name"]]
        for row in self.row_extractor():
            yield row[sample_col_idx]

    @staticmethod
    def get_file_type(file):
        p = subprocess.Popen(["file", "-b", file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        human, err = p.communicate()
        p = subprocess.Popen(["file","-b","-i", file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        mime, err = p.communicate()
        return human.decode().rstrip(), mime.decode().rstrip()

    @staticmethod
    def is_text_file(file):
        human, mime = PhenoReader.get_file_type(file)
        istext = re.search(r"\btext\b", human) is not None
        return istext, human, mime  

if __name__ == "__main__":
    import json
    import sys
    if len(sys.argv)>=2:
        meta = None
        if len(sys.argv)>=3:
            with open(sys.argv[2]) as f:
                meta = json.load(f)
            p = PhenoReader(sys.argv[1], meta)
            print([x for x in p.row_extractor()])
        else:
            p = PhenoReader(sys.argv[1], meta)
            print(json.dumps(p.meta, indent=2))

