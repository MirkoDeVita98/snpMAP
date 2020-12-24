import json 

def load_config():
    with open('config.json') as config_file:
        data = json.load(config_file)
        return data