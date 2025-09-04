import re
from itertools import chain

def sanitize(x):
    if re.match(r'^[^A-Za-z]', x):
        x = "X" + x
    x = re.sub(r'[^A-Za-z0-9._]', "_", x)
    return x

class ColumnFactory:
    @staticmethod
    def get_by_name(name, pr):
        if isinstance(name, dict):
            options = name
            name = name.get("name", "")
        else:
            options = None
        coldef = next((x for x in pr.meta["columns"] if x["name"]==name), None)
        return ColumnFactory.__get_column_class(coldef, pr, options) 

    @staticmethod
    def get_by_special_class(colclass, pr):
        coldef = next((x for x in pr.meta["columns"] if x["class"]==colclass), None)
        if colclass in ["family_id", "sample_id", "father_id", "mother_id", "sex"]:
            return PedRequiredColumn(coldef, colclass, pr)
        else:
            return None

    @staticmethod
    def find_id_column(cols):
        col = next((x for x in cols if x.coldef and x.coldef["class"] == "sample_id"), None)
        return col 

    @staticmethod
    def __get_column_class(coldef, pr, options=None):
        if not "class" in coldef:
            raise Exception("Invalid Column Definition")
        colclass = coldef["class"]
        if colclass=="categorical":
            return CategoricalColumn(coldef, pr, options=options)
        elif colclass=="binary":
            return BinaryColumn(coldef, pr, options=options)
        else:
            return Column(coldef, pr, options=options)

class Column(object):
    def __init__(self, coldef, pr, options=None):
        self.coldef = coldef
        self._raw_values = []
        if self.coldef:
            self.name = coldef["name"]
            self.colindex = [x["name"] for x in pr.meta["columns"]].index(self.name)
            self.missing = coldef.get("missing", None)
            self.colclass = coldef.get("class", None)

    def headers(self):
        return [self.name]

    def header_specs(self):
            # One header, numeric/default column
            return [{
                "raw": self.headers()[0],
                "meta": {"source": self.name, "class": self.colclass, "type": "numeric", "level": None}
            }]


    def append(self, row):
        val = row[self.colindex]
        if self.missing is not None and val==self.missing:
            val = None
        self._raw_values.append(val)
        return val

    def value(self, index):
        val = self._raw_values[index] 
        return val

    def values(self, index):
        val = self.value(index)
        if val is None:
            return None
        return [val]

    def reorder(self, indexlist):
        new_raw_values = []
        for idx in indexlist:
            if idx is not None:
                new_raw_values.append(self._raw_values[idx])
            else:
                new_raw_values.append(None)
        self._raw_values = new_raw_values

    def __len__(self):
        return len(self._raw_values)

class CategoricalColumn(Column):
    def __init__(self, coldef, pr, options=None):
        super(CategoricalColumn, self).__init__(coldef, pr, options)
        self.levels = coldef["levels"]
        self.ref_level = self.levels[0]
        self.contr_levels = self.levels[1:]

    def headers(self):
        return [self.name + "_" + x for x in self.contr_levels]


    def header_specs(self):
            # One-hot per contrast level (exclude reference)
            return [{
                "raw": f"{self.name}_{lvl}",
                "meta": {"source": self.name, "class": "categorical", "type": "categorical", "level": lvl}
            } for lvl in self.contr_levels]

    def values(self, index):
        val = self.value(index) 
        if val is None:
            return None
        ret = ["0"] * len(self.contr_levels)
        if val == self.ref_level:
            return ret
        for index, level in enumerate(self.contr_levels):
            if val == level:
                ret[index] = "1"
                return ret
        raise Exception("Found unexpected value ({}) in categorical column ()".format(val, self.name))

class BinaryColumn(Column):
    
    def __init__(self, coldef, pr, options=None):
        super(BinaryColumn, self).__init__(coldef, pr, options)
        self.levels = coldef["levels"]
        if options and "event" in options:
            self.set_event_level(options["event"])
        else:
            self.set_event_level(self.levels[-1])

    def headers(self):
        return [self.name + "_" + self.event_level]

    def header_specs(self):
           hdr = f"{self.name}_{self.event_level}"
           return [{
               "raw": hdr,
               "meta": {"source": self.name, "class": "binary", "type": "binary", "level": self.event_level}
           }]

    def values(self, index):
        val = self.value(index) 
        if val is None:
            return None
        if not val in self.levels:
            raise Exception("Found unexpected value in binary column")
        if val == self.event_level:
            return ["1"]
        else:
            return ["0"]

    def set_event_level(self, level):
        if not level in self.levels:
            raise Exception("Invalid event level: {}".format(level))
        self.event_level = level

class PedRequiredColumn(Column):
    def __init__(self, coldef, field, pr):
        super(PedRequiredColumn, self).__init__(coldef, pr)
        self.field = field
        print("filed is ",self.field)
        if not hasattr(self, "name") or self.name is None:
            self.name = field
        else:
            print("Using ped field {} mapped to column {}".format(field, self.name))

    def headers(self):
        header = ""
        if self.field == "family_id":
            header = "FAM_ID"
        elif self.field == "sample_id":
            header = "IND_ID"
        elif self.field == "father_id":
            header = "FAT_ID"
        elif self.field == "mother_id":
            header = "MOT_ID"
        elif self.field == "sex":
            header = "SEX"
        else:
            header = self.field
        return [header]


    def header_specs(self):
            ped_hdr = self.headers()[0]
            source = getattr(self, "name", None) or self.field
            return [{
                "raw": ped_hdr,
                "meta": {"source": source, "class": self.field, "type": "ped", "level": None}
            }]

    def append(self, row):
        if self.coldef is None or self.colindex is None:
            return None
        return super(PedRequiredColumn, self).append(row)

    def values(self, index):
        if self.coldef is None:
            return ["0"]
        return super(PedRequiredColumn, self).values(index)

def flatten(x):
    return list(chain.from_iterable((v if v is not None else [None] for v in x )))

class PedWriter:
    def __init__(self, phenoreader=None, resp=None, covar=None, samples=None):
        pedcols = ["family_id", "sample_id", "father_id", "mother_id", "sex"]
        self.pedcols = [ColumnFactory.get_by_special_class(x, phenoreader) for x in pedcols]
        self.respcols = [ColumnFactory.get_by_name(resp, phenoreader)]
        self.covarcols = [ColumnFactory.get_by_name(x, phenoreader) for x in covar]
        self.allcols = self.pedcols + self.respcols + self.covarcols
        for row in phenoreader.row_extractor(samples=samples):
            for col in self.allcols:
                col.append(row)
        self.expand_columns()

    def merge_covar(self, phenoreader=None, covar=None):
        datacols =  [ColumnFactory.get_by_name(x, phenoreader) for x in covar]
        matchcol = ColumnFactory.get_by_special_class("sample_id", phenoreader)
        if matchcol is None:
            raise Exception("Unable to find Sample ID column")
        lookup = dict()
        for idx, row in enumerate(phenoreader.row_extractor()):
            id = matchcol.append(row)
            lookup[id] = idx
            for col in datacols:
                col.append(row)

        reindex = []
        matchcol = ColumnFactory.find_id_column(self.allcols)
        if matchcol is None:
            raise Exception("Unable to find Sample ID column")
        for idx in range(len(matchcol)):
            sample = matchcol.value(idx)
            if sample in lookup:
                reindex.append(lookup[sample])
            else:
                reindex.append(None)

        for col in datacols:
            col.reorder(reindex)

        self.covarcols += datacols
        self.allcols += datacols
        self.expand_columns()

    def expand_columns(self):

      def uniqueify(vals, existing):
          taken = set(existing)
          uniqued = []
          for val in vals:
              newval = sanitize(val)
              ind = 1
              while newval in taken:
                  newval = sanitize(f"{val}.{ind}")
                  ind += 1
              uniqued.append(newval)
              taken.add(newval)
          assert len(vals) == len(uniqued)
          return uniqued

      def collect_specs(cols):
          specs = []
          for c in cols:
              specs.extend(c.header_specs())   # <--- new: each Column returns header + meta
          return specs

      ped_specs   = collect_specs(self.pedcols)
      resp_specs  = collect_specs(self.respcols)
      covar_specs = collect_specs(self.covarcols)

      self.header_meta = {}   # <--- new: dict of final header → metadata

      # Ped headers
      ped_raw = [s["raw"] for s in ped_specs]
      ped_uni = uniqueify(ped_raw, [])
      self.pedheaders = ped_uni
      for name, spec in zip(ped_uni, ped_specs):
          self.header_meta[name] = spec["meta"]

      self.headers = self.pedheaders[:]

      # Resp headers
      resp_raw = [s["raw"] for s in resp_specs]
      resp_uni = uniqueify(resp_raw, self.headers)
      self.respheaders = resp_uni
      self.headers += resp_uni
      for name, spec in zip(resp_uni, resp_specs):
          self.header_meta[name] = spec["meta"]

      # Covar headers
      print("Covar specs:", covar_specs)
      covar_raw = [s["raw"] for s in covar_specs]
      covar_uni = uniqueify(covar_raw, self.headers)
      print("Covar uni:", covar_uni)
      self.covarheaders = covar_uni
      self.headers += covar_uni
      for name, spec in zip(covar_uni, covar_specs):
          self.header_meta[name] = spec["meta"]

#     def expand_columns(self):
#
#         def uniqueify(vals, existing):
#             taken = existing[:]
#             uniqued = []
#             for val in vals:
#                 newval = sanitize(val)
#                 ind = 1
#                 while newval in taken:
#                     newval = sanitize(val + "." + ind)
#                 uniqued.append(newval)
#             assert len(vals) == len(uniqued)
#             return uniqued
#
#         self.pedheaders = uniqueify(flatten([x.headers() for x in self.pedcols]), [])
#         self.headers = self.pedheaders
#         self.respheaders = uniqueify(flatten([x.headers() for x in self.respcols]), self.headers)
#         self.headers = self.headers + self.respheaders
#         self.covarheaders = uniqueify(flatten([x.headers() for x in self.covarcols]), self.headers)
#         self.headers = self.headers + self.covarheaders

    def get_response_headers(self):
        return self.respheaders

    def get_covar_headers(self):
        return self.covarheaders
    def get_header_meta(self): return self.header_meta

    def get_covar_meta(self): return {h: self.header_meta[h] for h in self.covarheaders}

    def _result_rows(self):
        for idx in range(max((len(x) for x in self.allcols))):
            vals = flatten([x.values(idx) for x in self.allcols])
            has_missing = any((x is None for x in vals))
            if not has_missing:
                yield vals

    def get_response_values(self):
        id_index = 1
        resp_index = 5
        resp = dict()
        for vals in self._result_rows():
            resp[vals[id_index]] = float(vals[resp_index])
        return resp

    def write_to_file(self, fconn, comment_header=True):
        row_count = 0
        self.expand_columns()
        header = "\t".join(self.headers) + "\n"
        if comment_header:
            header = "#" + header
        fconn.write(header)
        for vals in self._result_rows():
            fconn.write("\t".join(vals) + "\n")
            row_count += 1
        return row_count

if __name__ == "__main__":
    from .pheno_reader import PhenoReader
    import json
    import os
    import sys

    class screenout:
        def write(self, txt):
            sys.stdout.write(txt)

    def init(filename):
        pr = PhenoReader(filename)
        meta = pr.infer_meta()
        print(json.dumps(meta))

    def get_pr(filename):
        meta = None
        with open(filename.replace(".txt",".json")) as f:
            meta = json.load(f)
        return PhenoReader(filename, meta)

    #init(os.path.expanduser("~/in2.txt"))
    pr1 = get_pr(os.path.expanduser("~/in1.txt"))
    pr2 = get_pr(os.path.expanduser("~/in2.txt"))
    pedw = PedWriter(pr1, "col3",["col4","col2"])
    sout = screenout();
    pedw.write_to_file(sout)
    pedw.merge_covar(pr2, ["pc3", "pc1"])
    pedw.write_to_file(sout)
