import os
import configparser

config = configparser.ConfigParser()

config.read('somatic_config.ini')

print(config.sections())

config['test'] = {
    'a':'b',
    'path':os.environ['PWD'],
}
#config['test']['a'] = 'b'

with open('new_config.ini','w') as odata:
    config.write(odata)