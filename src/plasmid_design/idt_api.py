from base64 import b64encode
import json
from urllib import request, parse
import requests
import joblib
from .utils import ProgressParallel


# IDT API reference
# https://www.idtdna.com/pages/tools/apidoc
# https://www.idtdna.com/Restapi/swagger/ui/index#/Complexity


idt_username = 'bilforama'
idt_password = 'd256siaGBDze$jC'
client_id = 'dsl389ds'
client_secret = '920d3cce-5580-47f9-b4ed-bc5d5837d2bc'

max_parallel_idt_queries = 50

# scores for gblocks and genes can be quite different
idt_api_endpoints = {
    'gblock_complexity': 'https://www.idtdna.com/Restapi/v1/Complexities/ScreenGblockSequences',
    'gene_complexity': 'https://www.idtdna.com/Restapi/v1/Complexities/ScreenGeneSequences',
}


def check_dna_idt_parallel(list_of_dna, method, progress=True):
    """
    """

    token = get_access_token()
    def f(dna):
        response = check_dna_idt([dna], token=token, method=method, text_only=False)
        seq = json.loads(response.request.body)[0]['Sequence']
        return (seq, json.loads(response.text))

    parallel = ProgressParallel if progress else joblib.Parallel

    issues = dict(parallel(n_jobs=max_parallel_idt_queries, prefer='threads')(
        [joblib.delayed(f)(x) for x in list_of_dna]
    ))
    return [issues[x][0] for x in list_of_dna]


def score_dna_idt_parallel(list_of_dna, method):
    results = check_dna_idt_parallel(list_of_dna, method)
    return [extract_score(x) for x in results]


def extract_score(issues):
    return sum(x['Score'] for x in issues)    


def check_dna_idt(list_of_dna, method, token=None, text_only=True):
    """Results are sometimes out-of-order if more than one sequence is tested??
    :param genes: list of DNA
    :param method: either "gblock" or "gene"
    """
    if 1 < len(list_of_dna):
        raise NotImplementedError('Request one sequence at a time or '
                                  'use check_dna_idt_parallel')

    data = format_genes_for_idt_api(enumerate(list_of_dna))
    if token is None:
        token = get_access_token()


    url = idt_api_endpoints[method + '_complexity']

    headers = {
        'Content-Type': 'application/json', 
        'Authorization': f'Bearer {token}',
    }
    response = requests.request("POST", url, headers=headers, data=data)
    if text_only:
        response = json.loads(response.text)
    return response


def format_genes_for_idt_api(genes):
    """Create message body
    """
    xs = [{'Name': a, 'Sequence': b} for a,b in genes]
    return str(xs).replace("'", '"')


def get_access_token():
    return get_access_token_idt(client_id, client_secret, idt_username, idt_password)


def get_access_token_idt(client_id, client_secret, idt_username, idt_password):
    """
    From IDT

    Create the HTTP request, transmit it, and then parse the response for the 
    access token.
    
    The body_dict will also contain the fields "expires_in" that provides the 
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]


