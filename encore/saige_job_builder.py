
from .base_model import BaseModel
from .ped_writer import PedWriter
import os
from .chunk_progress import get_chr_chunk_progress
import yaml


class SaigeModel(BaseModel):
    filters = [("min-mac-20", "MAC > 20"),
        ("min-maf-001", "MAF > 0.1%"),
        ("min-maf-001-mac-20","MAF > 0.1% AND MAC > 20")]

    def __init__(self, working_directory="./", app_config=None):
        BaseModel.__init__(self, working_directory, app_config) 
        self.cores_per_job = 56

    def get_opts(self, model, geno):
        opts = {}
        if model.get("response_invnorm", False):
             opts['inv_norm']= 'true'
        if model.get("variant_filter", False):
            vf = model.get("variant_filter")
            if vf == "min-maf-001-mac-20":
                opts['min_mac']=20
                opts['min_maf']= 0.001
            elif vf == "min-maf-001":
                #opts.append("STEP2OPT='--minMAF 0.001 --IsOutputAFinCaseCtrl=FALSE'")
                opts['min_maf']= 0.001
            elif vf == "min-mac-20":
                opts['min_mac']= 20

            else:
                raise Exception("Unrecognized variant filter ({})".format(vf))

        opts['region_size']=10000000
        if model.get("region", None):

                    region = model.get("region").upper()
                    print("reiong", region)
                    if region.startswith("CHR"):
                        region = region[3:]
                    opts["CHRS"] = region

        # elif geno.get_chromosomes():
        #     opts.append("CHRS='{}'".format(geno.get_chromosomes()))
        return opts

    def get_analysis_commands(self, model_spec, geno, ped,pheno):


        pipeline = self.app_config["SAIGE_SIF_FILE"][0]

        if "SAIGE_BINARY" in self.app_config:
            binary = self.app_config["SAIGE_BINARY"]
        if isinstance(binary, tuple):
            binary = binary[0]
        if not binary:
            raise Exception("Unable to find Saige sif file file  (pipeline: {})".format(pipeline))
        #for dev singularity exec -B /net/wonderland:/net/wonderland:ro,/net/dumbo:/net/dumbo:ro
        #cmd = "singularity exec -B /net/encore1/savant:/net/encore1/savant:ro -B /net/encore1/encoredata:/net/encore1/encoredata {} ".format(pipeline) + \
        cmd = "singularity exec -B /net/wonderland:/net/wonderland:ro  -B /net/dumbo:/net/dumbo {} ".format(pipeline) + \
              " snakemake --snakefile {}".format(binary)+ \
              " -j ${SLURM_CPUS_PER_TASK}"

        optlist = self.get_opts(model_spec, geno)
        optlist["post_processing_script_dir"]=self.app_config["POST_PROCESSING_SCRIPT"]
        optlist["gene_annotation_bed"]=self.app_config["NEAREST_GENE_BED"]
        optlist["savs_path"]= geno.get_sav_path(1).replace("chr1.sav", "")
        optlist["plinkFile"]= geno.get_pca_genotypes_path()
        optlist["samples_file"]= self.working_directory + "/samples.txt"

        optlist["output_dir"]=self.working_directory
        optlist["phenoFile"]= ped.get("path")
        optlist["outputPrefix"]= "step1"
        optlist["nThreads"]= 54
        optlist["sampleIDColinphenoFile"]= "IND_ID"
        optlist["IsOverwriteVarianceRatioFile"]= "TRUE"
        print("just above **********",self.app_config["BIND_MOUNT"])
        optlist["BIND_MOUNT"]=self.app_config["BIND_MOUNT"]


        resplist = ped.get("response")

        if len(resplist)>0:
                optlist['phenoCol']=resplist[0]
        covars = ped.get("covars")
        covar_meta = ped.get("covar_meta") or {}
        if len(covars)>0:
            #cmd += " COVAR={}".format(",".join(covars))
            print("covars are",covars)
            cat_covars = [h for h in covars
                          if covar_meta.get(h, {}).get("type") in {"categorical","binary"}]
            print("cat_covars",cat_covars)
            if covars:
                optlist['covarColList']= "{}".format(",".join(covars))          # → map to SAIGE --covarColList
            if cat_covars:
                optlist['qCovarColList']="{}".format(",".join(cat_covars))



        confilepath = self.relative_path("config.yml")
        with open(confilepath, 'w') as file:
            documents = yaml.dump(optlist, file)
        return [cmd]

    def get_postprocessing_commands(self, geno, result_file="./results.txt.gz"):
        cmds = []
        cmds.append("zcat -f {} | ".format(result_file) + \
            'awk -F"\\t" \'BEGIN {OFS="\\t"} NR==1 {for (i=1; i<=NF; ++i) {if($i=="p.value") pcol=i; if($i=="N") ncol=i}; if (pcol<1 || ncol<1) exit 1; print} ' + \
            '($ncol > 0 && $pcol < 0.001) {print}\' | ' + \
            "{} -c > output.filtered.001.gz".format(self.app_config.get("BGZIP_BINARY", "bgzip")))
        if self.app_config.get("MANHATTAN_BINARY"):
            cmd  = "{} {} ./manhattan.json".format(self.app_config.get("MANHATTAN_BINARY", ""), result_file)
            cmds.append(cmd)
        if self.app_config.get("TOPHITS_BINARY"):
            cmd = "{} {} ./qq.json".format(self.app_config.get("QQPLOT_BINARY", ""), result_file)
            cmds.append(cmd)
        if self.app_config.get("TOPHITS_BINARY"):
            cmd =  "{} {} ./tophits.json".format(self.app_config.get("TOPHITS_BINARY"), result_file)
            if geno.get_build_nearest_gene_path():
                cmd += " --gene {}".format(geno.get_build_nearest_gene_path())
            cmds.append(cmd)
        return cmds 

    def write_ped_file(self, ped_file_path, model_spec, geno=None, pheno=None):
        geno = self.get_geno(model_spec, geno)
        pheno = self.get_pheno(model_spec, pheno)
       
        try:
            ped_writer = self.get_ped_writer(model_spec, geno, pheno) 
            with open(ped_file_path, "w") as pedfile:
                ped_writer.write_to_file(pedfile, comment_header=False)
            return {
                "response": ped_writer.get_response_headers(),
                "covars": ped_writer.get_covar_headers(),
                "covar_meta": getattr(ped_writer, "get_covar_meta", lambda: {})(),  # ✅ add this
                "path": ped_file_path
            }
        except Exception as e:
            raise Exception("Failed to create ped file ({})".format(e))



    def prepare_job(self, model_spec):
        geno = self.get_geno(model_spec)
        pheno = self.get_pheno(model_spec)

        ped = self.write_ped_file(self.relative_path("pheno.ped"), model_spec, geno, pheno)
        cmds =  self.if_exit_success(
            self.get_analysis_commands(model_spec, geno, ped,pheno),
            self.get_postprocessing_commands(geno))
        return {"commands": cmds}

    def get_progress(self):
        output_file_glob = self.relative_path("step2.bin.*.txt")
        fre = r'step2\.bin\.(?P<chr>\w+)\.txt$'
        resp = get_chr_chunk_progress(output_file_glob, fre)
        return resp

    def validate_model_spec(self, model_spec):
        if not "pipeline_version" in model_spec:
            if "SAIGE_VERSION" in self.app_config:
                model_spec["pipeline_version"] = self.app_config["SAIGE_VERSION"]
        if "variant_filter" in model_spec:
            if hasattr(self, "filters"):
                vf = model_spec["variant_filter"]
                if not any([vf == x[0] for x in self.filters]):
                    raise Exception("Did not recognize variant filter ({})".format(vf))
            else:
                raise Exception("No variant filters defined but one was requested ({})".format(e))
        else:
            if hasattr(self, "filters"):
                model_spec["variant_filter"] = self.filters[0][0]

    def get_failure_reason(self):
        log_file_path = self.relative_path("saige.log")
        if not os.path.isfile(log_file_path):
            return None
        with open(log_file_path, 'rt') as f:
            for line in f:
                if "matrix is singular" in line:
                    return "Matrix is singular or not positive definite"
                elif "parameter estimate is 0" in line:
                                    return "The first variance component parameter estimate is 0"
                elif "variance of the phenotype is much smaller" in line:
                    return ("Variance of the phenotype is much smaller than 1. "
                            "Please consider using inverse normalized response")
        return None
        
class LinearSaigeModel(SaigeModel):
    model_code = "saige-qt"
    model_name = "Saige Linear Mixed Model"
    model_desc = "Fast linear mixed model with kinship adjustment"
    depends = ["sav"]

    def __init__(self, working_directory, config):
        SaigeModel.__init__(self, working_directory, config) 

    def get_opts(self, model, geno):
        opts = super(self.__class__, self).get_opts(model, geno) 
        opts["traitType"]="quantitative"
        opts["model"]=self.model_code
        return opts

class BinarySaigeModel(SaigeModel):
    model_code = "saige-bin"
    model_name = "Saige Logistic Mixed Model"
    model_desc = "Fast logistic regression model with kinship adjustment"
    depends = ["sav", "snps"]
    response_class = "binary"

    def __init__(self, working_directory, config):
        SaigeModel.__init__(self, working_directory, config)

    def get_opts(self, model, geno):
        opts = super(self.__class__, self).get_opts(model, geno) 
        opts["traitType"]="binary"
        opts["model"]=self.model_code
        return opts

    def validate_model_spec(self, model_spec):
        super(self.__class__, self).validate_model_spec(model_spec) 
        resp = model_spec.get("response", None)
        if not resp:
            raise Exception("Missing model response")

        if not isinstance(resp, dict):
            resp = {"name": resp}
        resp_name = resp.get("name", None)
        resp_event = resp.get("event", None)

        pheno = self.get_pheno(model_spec)
        levels = pheno.get_column_levels(resp_name)
        if not levels:
            raise Exception(("Variable {} does not appear to be a binary trait "
                "(No levels found)").format(resp_name))
        if len(levels) != 2:
            raise Exception("Response must have exactly 2 levels, found {}: {}".
                format(len(levels), ",".join(levels)))
        if resp_event:
            if not resp_event in levels:
                raise Exception(("Requested event level does not match data. "
                    "Requested: {}; Data: {}").format(resp_event, ",".join(levels)))
        else:
            resp_event = levels[1]

        resp["event"] = resp_event
        model_spec["response"] = resp
        model_spec["response_desc"] = "Pr({} = {})".format(resp_name, resp_event)

