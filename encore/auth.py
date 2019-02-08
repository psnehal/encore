import inspect
from functools import wraps
from .job import Job
from .phenotype import Phenotype
from flask import current_app
from flask_login import current_user

#most of these are helper function
#there are the dectorators you should use
# - check_view_job (needs job_id and job=None parameters)
#      will inject job parameter if None
#      checks if current user has access to job data
# - check_edit_job (needs job_id and job=None parameters)
#      will inject job parameter if None
#      checks if current user can update job data
# - access_pheno_page (needs pheno_id and pheno=None parameters)
#      will inject pheno parameter if None
#      checks if current user has access to phenotype data
# - check_edit_pheno (needs pheno_id and pheno=None)

def can_user_view_job(user, job):
    if not user or not job:
        return False
    user_id = user.rid
    if user_id == job.user_id:
        return True
    if user_id in (x["user_id"] for x in job.users):
        return True
    if user.is_admin():
        return True
    return False

def can_user_edit_job(user, job):
    if not user or not job:
        return False
    user_id = user.rid
    if user_id == job.user_id:
        return True
    if user.is_admin():
        return True
    return False

def can_user_view_pheno(user, pheno):
    if not user or not pheno:
        return False
    user_id = user.rid
    if user_id == pheno.user_id:
        return True
    if user.is_admin():
        return True
    return False

def can_user_edit_pheno(user, pheno):
    if not user or not pheno:
        return False
    user_id = user.rid
    if user_id == pheno.user_id:
        return True
    if user.is_admin():
        return True
    return False

# decorators...
def splat_args(f, asf):
    formals, pargs, pkwargs, defaults = inspect.getargspec(asf)
    @wraps(f)
    def inner(*args, **kwargs):
        kwargs.update({k:v for (k,v) in zip(formals, args)})
        return f(**kwargs)
    return inner

def inject_job(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        if "job_id" in kwargs and not kwargs.get("job",None):
            kwargs["job"] = Job.get(kwargs["job_id"], current_app.config)
        else:
            raise Exception("required 'job_id' parameter not found")
        return f(**kwargs)
    return inner

def can_view_job(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        job = kwargs["job"]
        user = current_user
        if job is not None:
            if can_user_view_job(user, job):
                return f(**kwargs)
            else:
                return "Unauthorized", 403
        else:
            return "Job not found", 404
    return inner

def check_view_job(f):
    return splat_args(inject_job(can_view_job(f)), f)

def can_edit_job(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        job = kwargs["job"]
        user = current_user
        if job is not None:
            if can_user_edit_job(user, job):
                return f(**kwargs)
            else:
                return "Unauthorized", 403
        else:
            return "Job not found", 404
    return inner

def check_edit_job(f):
    return splat_args(inject_job(can_edit_job(f)), f)

def inject_pheno(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        if "pheno_id" in kwargs and not kwargs.get("pheno",None):
            kwargs["pheno"] = Phenotype.get(kwargs["pheno_id"], current_app.config)
        else:
            raise Exception("required 'pheno_id' parameter not found")
        return f(**kwargs)
    return inner

def can_view_pheno_page(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        pheno = kwargs["pheno"]
        user = current_user
        if pheno is not None:
            if can_user_view_pheno(user, pheno):
                return f(**kwargs)
            else:
                return "Unauthorized", 403
        else:
            return "Phenotype not found", 404
    return inner

def access_pheno_page(f):
    return splat_args(inject_pheno(can_view_pheno_page(f)), f)

def can_edit_pheno(f):
    @wraps(f)
    def inner(*args, **kwargs):
        if len(args):
            raise Exception("unnamed args")
        pheno = kwargs["pheno"]
        user = current_user
        if pheno is not None:
            if can_user_edit_pheno(user, pheno):
                return f(**kwargs)
            else:
                return "Unauthorized", 403
        else:
            return "Phenotype not found", 404
    return inner

def check_edit_pheno(f):
    return splat_args(inject_pheno(can_edit_pheno(f)), f)

def admin_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if not current_user.is_admin():
            return "You do not have access", 403 
        return f(*args, **kwargs)
    return decorated_function
