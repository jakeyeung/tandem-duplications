'''
Created on 2013-05-17

@author: jyeung
'''


import sys
import time
import os


_cur_dir = os.path.dirname(os.path.realpath(__file__))
_upthree_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
_plot_dir = os.path.join(_upthree_dir, 'inputs')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Tandem dupes file and TSS data must be given on the command line.')
        sys.exit()
    starttime = time.time()
    tandem_dupe_fname = sys.argv[1]
    tss_data_fname = sys.argv[2]
    print _cur_dir, _upthree_dir, _plot_dir
    print(time.time() - starttime)