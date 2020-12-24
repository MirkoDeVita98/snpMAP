import os.path
import urllib.request
import zipfile

from utils import parser

config = parser.load_config()
cache_dir = config['output']['cache_dir']


# Most likely not needed
def download_SNPEff():
    snpeff_url = 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip'
    download_fpath = os.path.join(cache_dir, 'snpEff_latest_core.zip')
    if not os.path.exists(download_fpath):
        urllib.request.urlretrieve(snpeff_url, download_fpath)
    
    if not os.path.exists('dependencies/snpEff'):
        with zipfile.ZipFile(download_fpath, 'r') as zip_ref:
            zip_ref.extractall(cache_dir)