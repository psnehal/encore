from flask import Blueprint, render_template, request, json, current_app, redirect, session, url_for, redirect
from flask_login import LoginManager, logout_user, login_required, current_user
import requests
from rauth import OAuth2Service
from .user import User
import flask_login
from . import sql_pool
import jwt
import datetime
from pprint import pprint
from urllib.parse import urlencode
from oic.oic.message import RegistrationResponse
from oic.utils.authn.client import CLIENT_AUTHN_METHOD
from oic.oic.message import ProviderConfigurationResponse



# or
# op_info = ProviderConfigurationResponse(**info)
# if you have the provider info in the form of a dictionary




import os
import requests

from oic import rndstr
from oic.oauth2 import AuthorizationResponse
from oic.oic import Client
from oic.oic.message import OpenIDSchema, AccessTokenResponse
from oic.utils.authn.client import ClientSecretBasic, ClientSecretPost

googleinfo = requests.get("https://accounts.google.com/.well-known/openid-configuration")
google_params = googleinfo.json()


umichinfo=requests.get("https://shibboleth.umich.edu/.well-known/openid-configuration")
#umichinfo=requests.get("https://shib-idp-test.www.umich.edu/.well-known/openid-configuration")

umich_params=umichinfo.json()

# AUTHMACHINE_URL = "https://shib-idp-test.www.umich.edu/.well-known/openid-configuration"
# AUTHMACHINE_API_TOKEN = umich_params.get("token_endpoint")
ENCORE_SCOPE = 'openid edumember email'


auth = Blueprint("auth", __name__)

login_manager = LoginManager()

def __init__(self):
    self.client = self.get_client()
    if request.is_secure:
        proto = 'https://'
    else:
        proto = 'http://'
    self.host = proto + request.host
    print(self.host)

def encode_auth_token(user_id):
    try:
        payload = {
            'exp': datetime.datetime.utcnow() + datetime.timedelta(days=7, seconds=0),
            'iat': datetime.datetime.utcnow(),
            'sub': user_id
        }
        return jwt.encode(
            payload,
            current_app.config.get('JWT_SECRET_KEY'),
            algorithm='HS256'
        )
    except Exception as e:
        print(e)
        return None

def decode_auth_token(auth_token):
    try:
        payload = jwt.decode(auth_token, current_app.config.get('JWT_SECRET_KEY'))
        return payload['sub']
    except jwt.ExpiredSignatureError:
        return 'Signature expired. Please log in again.'
    except jwt.InvalidTokenError:
        return 'Invalid token. Please log in again.'

@login_manager.request_loader
def user_loader_from_request(request):
    auth_header = request.headers.get("Authorization")
    if auth_header and auth_header.startswith("Bearer "):
        auth_token = auth_header.split(" ")[1]
    else:
        auth_token = ""
    if auth_token:
        try:
            email = decode_auth_token(auth_token)
            return load_user(email)
        except Exception as e:
            print(e)
    return None

@login_manager.user_loader
def user_loader(email):
    return load_user(email)

@auth.route("/sign-in", methods=["GET"])
def get_sign_in():
    if request.args.get("orig", None):
        session["post_login_page"] = request.args.get("orig")
    return get_sign_in_view("sign-in")

@auth.route("/signin", methods=["GET"])
def get_signin():
     if request.args.get("orig", None):
         session["post_login_page"] = request.args.get("orig")
     return get_check_in_oidcview("signin")
    #return redirect(get_authorization_url())


@auth.route("/checkin", methods=["GET"])
def get_checkin():
    if request.args.get("orig", None):
        session["post_login_page"] = request.args.get("orig")
    return get_check_in_oidcview("checkin")


def get_check_in_oidcview(target):
    #print("sign is oidc : "+target)
    signin_url = request.url_root + target
    print("sign is url : "+signin_url)
    client = Client(client_authn_method=CLIENT_AUTHN_METHOD)
    info = {"client_id": current_app.config.get("REMOVED", None), "client_secret": current_app.config.get("REMOVED", None)}
    client_reg = RegistrationResponse(**info)
    client.store_registration_info(client_reg)

    # #print(request.url_root + 'oidc-callback')
    # print(umich_params.get("authorization_endpoint"))
    # print(umich_params.get("token_endpoint"))
    # print(umich_params.get("userinfo_endpoint"))
    # print(umich_params.get("jwks_uri"))

    op_info = ProviderConfigurationResponse(
        version="1.0", issuer=umich_params.get("issuer"),
        authorization_endpoint=umich_params.get("authorization_endpoint"),
        token_endpoint=umich_params.get("token_endpoint"),
        userinfo_endpoint=umich_params.get("userinfo_endpoint"),
        jwks_uri=umich_params.get("jwks_uri"),
        )
    client.handle_provider_config(op_info, op_info['issuer'])

    # //SO if the user is part of that encore group, then will check if the user exist in the table using unique name , if not then create user using following command
    # //and then allow acccess
    if "code" in request.args:
        authorization_response = client.parse_response(
            AuthorizationResponse,
            info=request.args,
            sformat='dict')
        print("*****************************************************authorization reposnse***************************")
        args = {
            'code': authorization_response['code'],
            'client_id': client.client_id,
            'client_secret': client.client_secret,
            'redirect_uri': signin_url
        }

        access_token_ret= client.do_access_token_request(
            scope=ENCORE_SCOPE,
            state=authorization_response['state'],
            request_args=args,
            authn_method='client_secret_post')

        print("***************************************************** access_token_ret ***************************")
        print(access_token_ret)
        access_token2 = access_token_ret['access_token']

        #userinfo_request(access_token)
        # Parameters:	access_token (str) – Bearer access token to use when fetching userinfo
        # Returns:	UserInfo Response
        # Return type:	oic.oic.message.OpenIDSchema

        user_info2=client.user_info_request(
            access_token=access_token2
        )
        print("***************************************************** user_info2 ***************************")
        print(user_info2)

        print("***************************************************** Running user_info now  ***************************")
        user_info = client.do_user_info_request(
            access_token=access_token2)
        print("***************************************************** user_info ***************************")
        print(user_info)
        print(user_info['sub'])
        useremail= user_info['email']
        usersub = user_info['sub']
        groupinfo= user_info['edumember_ismemberof']
        user = load_uniquename(usersub)

        if 'encore mgi' in groupinfo:
            print("user is inside the group")
            #print("***************** is present or not in the db")
            if user:
                print("user is present in the db")
                if user.is_active():
                    print("user is active")
                    flask_login.login_user(user)
                    redirect_to = session.pop("post_login_page", None)
                    print(redirect_to)
                    try:
                        endpoint, arguments = current_app.url_map.bind('localhost').match(redirect_to)
                        #print("in try")
                        print(endpoint)
                        print(arguments)
                    except Exception as e:
                        redirect_to = None
                    if redirect_to:
                        return redirect(redirect_to)
                    else:
                        #print("user is not active")
                        return redirect(url_for("user.index"))
                else:
                    error_message = "Account not active ({})".format(user_info["email"])
                    return render_template("/signin.html", error_message=error_message)
            else:
                db = sql_pool.get_conn()
                print("user is not present in the db")
                print("create user")
                userdev= {}
                userdev['email']=useremail
                userdev['fullname']=usersub
                userdev['uniquename']=usersub
                userdev['affiliation']='test'
                #userdev['creation_date']="DATE_FORMAT(2020-02-17 19:45:17, '%%Y-%%m-%%d %%H:%%i:%%s')"
                #userdev['last_login_date']="DATE_FORMAT(2020-02-17 19:45:17, '%%Y-%%m-%%d %%H:%%i:%%s')"
                userdev['can_analyze']=1
                userdev['is_active']=1

                ##Pass the dictionary to db cursor and then create user
                usercreate = User.createUser(userdev,db)
                usercreate.log_login(db)
                flask_login.login_user(usercreate)
                redirect_to = session.pop("post_login_page", None)

                try:
                    endpoint, arguments = current_app.url_map.bind('localhost').match(redirect_to)
                except Exception as e:
                    redirect_to = None
                if redirect_to:
                    return redirect(redirect_to)
                else:
                    return redirect(url_for("user.index"))
        else:
            print("user is not the part of the group")
            if user:
                print("user is not in the encore group but available in the database")
                user.is_active = 0
                db = sql_pool.get_conn()
                db.commit()
                error_message = "User is not part of the Encore_mgi mcommunity group ({})".format(user_info['email'])
                return render_template("/signin.html", error_message=error_message)
            error_message = "Not an authorized user ({})".format(user_info['email'])
            return render_template("/signin.html", error_message=error_message)
    elif "authorize" in request.args:
        print("authorize")
        nonce = rndstr()
        args = {
            'client_id': current_app.config.get("REMOVED", None),
            'response_type': 'code',
            'scope': ENCORE_SCOPE,
            'nonce': nonce,
            'redirect_uri': signin_url,
            'state': 'some-state-which-will-be-returned-unmodified'
        }
        url = umich_params.get("authorization_endpoint") + '?' + urlencode(args, True)
        #print(url)
        return redirect(url)
    else:
        return render_template("/signin.html")


@auth.route("/get-auth-token", methods=["GET"])
@login_required
def get_auth_token():
    return encode_auth_token(current_user.email)

@login_manager.unauthorized_handler
def unauthorized():
    if request.path.startswith("/api"):
        return "UNAUTHORIZED", 401
    else:
        orig = request.full_path
        if orig == "/?":
            orig = None
        return redirect(url_for("auth.get_signin", orig=orig))

@auth.route("/sign-out", methods=["GET"])
def sign_out():
    logout_user()
    return redirect(url_for("auth.get_signin"))



def load_userfullname(fullname):
    print("in load user fullname")
    print(fullname)
    db = sql_pool.get_conn()
    print(db)
    user = User.from_full_name(fullname, db)
    if user:
        #try:
        user.log_login(db)
        #except:
        #    pass
        return user
    else:
        return None



def load_uniquename(uniquename):
    print("in load user uniquename")
    print(uniquename)
    db = sql_pool.get_conn()
    print(db)
    user = User.from_unique_name(uniquename, db)
    if user:
        #try:
        user.log_login(db)
        #except:
        #    pass
        return user
    else:
        return None

def load_user(email):
    db = sql_pool.get_conn()
    user = User.from_email(email, db)
    if user:
        #try:
        user.log_login(db)
        #except:
        #    pass
        return user
    else:
        return None
def get_sign_in_oidcview(target):
    print("sign is oidc"+target)
    signin_url = request.url_root + target
    print("sign is url"+signin_url)

    # //SO if the user is part of that encore group, then will check if the user exist in the table using unique name , if not then create user using following command
    # //and then allow acccess

    db = sql_pool.get_conn()

    user = User.from_email('inspiresnehal@gmail.com', db)
    if user:
        print("found user")
    else:
        db = sql_pool.get_conn()
        print("user is not present in the db")
        print("create user")
        userdev= {}
        userdev['email']=""
        userdev['fullname']=""
        userdev['uniquename']='test'
        userdev['affiliation']=''
        userdev['creation_date']="DATE_FORMAT(jobs.creation_date, '%%Y-%%m-%%d %%H:%%i:%%s'), '%%Y-%%m-%%d %%H:%%i:%%s')"
        userdev['last_login_date']="DATE_FORMAT(2020-02-17 19:45:17, '%%Y-%%m-%%d %%H:%%i:%%s')"
        userdev['can_analyze']=0
        userdev['is_active']=0
        print("After the creating dictonary")
        usercreate = User.createUser(userdev,db)
        print("After the creating dictonary")

    oauth_service = OAuth2Service(
        name="Encore",
        client_id=current_app.config.get("REMOVED", None),
        client_secret=current_app.config.get("REMOVED", None),
        authorize_url=umich_params.get("authorization_endpoint"),
        base_url=umich_params.get("userinfo_endpoint"),
        access_token_url=umich_params.get("token_endpoint"))
    if "code" in request.args:
        print("1")
        print(request.args)
        oauth_session = oauth_service.get_auth_session(
            data={"code": request.args["code"],
                  "grant_type": "authorization_code",
                  "redirect_uri": signin_url},
            decoder=json.loads)

        # token = oauth_session.fetch_token(token_url=umich_params.get("token_endpoint"),
        #                           client_id=current_app.config.get("REMOVED", None),
        #                           client_secret=current_app.config.get("REMOVED", None),
        #                           include_client_id=True)

        print(dir(oauth_session))
        pprint(vars(oauth_session))
        print("access token")
        print(oauth_session.get("access_token"))

        print(oauth_session)
        r = oauth_session.get(umich_params.get("userinfo_endpoint"), params={'format': 'json'})
        print (r.json())
        for i, tweet in enumerate(r.json(), 1):
            print(tweet)
        #user_data = oauth_session.get("").json()
        #print(user_data)

        user = load_userfullname('snehal')
        if user:
            if user.is_active():
                print("user is active")
                flask_login.login_user(user)
                redirect_to = session.pop("post_login_page", None)
                try:
                    endpoint, arguments = current_app.url_map.bind('localhost').match(redirect_to)
                    print("in try")
                    print(endpoint)
                    print(arguments)
                except Exception as e:
                    redirect_to = None
                if redirect_to:
                    return redirect(redirect_to)
                else:
                    return redirect(url_for("user.index"))
            else:
                error_message = "Account not active ({})".format(user_data["email"])
                return render_template("/signin.html", error_message=error_message)
        else:
            error_message = "Not an authorized user ({})".format(user_data["email"])
            return render_template("/signin.html", error_message=error_message)
    elif "authorize" in request.args:
        print("authorize")
        authorize_url =oauth_service.get_authorize_url(
            scope=ENCORE_SCOPE,
            response_type="code",
            prompt="select_account",
            redirect_uri=signin_url)
        print(authorize_url)
        return redirect(authorize_url)
    else:
        return render_template("/signin.html")


def get_sign_in_view(target):
    signin_url = request.url_root + target
    oauth_service = OAuth2Service(
        name="google",
        client_id=current_app.config.get("GOOGLE_LOGIN_CLIENT_ID", None),
        client_secret=current_app.config.get("GOOGLE_LOGIN_CLIENT_SECRET", None),
        authorize_url=google_params.get("authorization_endpoint"),
        base_url=google_params.get("userinfo_endpoint"),
        access_token_url=google_params.get("token_endpoint"))

    if "code" in request.args:
        oauth_session = oauth_service.get_auth_session(
            data={"code": request.args["code"],
                  "grant_type": "authorization_code",
                  "redirect_uri": signin_url},
            decoder=json.loads)


        user_data = oauth_session.get("").json()
        print(user_data)
        user = load_user(user_data["email"])
        print(user)
        if user:
            if user.is_active():
                print("user is active")
                flask_login.login_user(user)
                redirect_to = session.pop("post_login_page", None)
                try:
                    endpoint, arguments = current_app.url_map.bind('localhost').match(redirect_to)
                except Exception as e:
                    redirect_to = None
                if redirect_to:
                    return redirect(redirect_to)
                else:
                    return redirect(url_for("user.index"))
            else:
                error_message = "Account not active ({})".format(user_data["email"])
                return render_template("/sign_in.html", error_message=error_message)
        else:
            error_message = "Not an authorized user ({})".format(user_data["email"])
            return render_template("/sign_in.html", error_message=error_message)
    elif "authorize" in request.args:
        authorize_url = oauth_service.get_authorize_url(
            scope="email",
            response_type="code",
            prompt="select_account",
            redirect_uri=signin_url)
        print(authorize_url)

        return redirect()
    else:
        return render_template("/sign_in.html")

