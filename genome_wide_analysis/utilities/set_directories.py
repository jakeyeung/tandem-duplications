'''
Created on 2013-05-19

@author: jyeung

Set directories, used at beginning of each main script. 
'''

import os


def set_directories(input_folder_name, output_folder_name, plot_folder_name):
    '''
    Return current directory of MAIN (not utilities)
    Get project root dir
    Get input dir
    Get output dir
    Get plot dir
    '''
    _cur_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    _proj_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
    _input_dir = os.path.join(_proj_dir, input_folder_name)
    _output_dir = os.path.join(_proj_dir, output_folder_name)
    _plot_dir = os.path.join(_output_dir, plot_folder_name)
    
    return _cur_dir, _proj_dir, _input_dir, _output_dir, _plot_dir

def set_input_dir(input_folder_name):
    _cur_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    _proj_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
    _input_dir = os.path.join(_proj_dir, input_folder_name)
    
    return _input_dir