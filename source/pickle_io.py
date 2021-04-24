import pickle
import os

def save_obj(obj, obj_name, variable_name):
    if not os.path.isdir(obj_name):
        os.mkdir(obj_name)
    v_path = obj_name + '/'+ variable_name + '.pkl'
    if os.path.exists(v_path):
        print('Warning: variable ' + variable_name + ' exists, overwritten')
    with open(v_path, 'wb') as f: # open file for writing in binary mode, delete if exists
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(obj_name, variable_name, **kwargs):
    v_path = obj_name + '/'+ variable_name + '.pkl'
    with open(v_path, 'rb') as f:
        if 'encoding' in kwargs:
            encoding = kwargs.get('encoding')
        else:
            encoding = 'ascii'
        return pickle.load(f, encoding = encoding)
