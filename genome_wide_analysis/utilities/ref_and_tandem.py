'''
Created on 2013-05-25

@author: jyeung
'''


import csv

        
class data(object):
    '''
    Reads from genome info file and tandem file.  
    '''
    
    def __init__(self, ref_path, tandem_path):
        '''
        Filepath
        '''
        self.refpath = ref_path
        self.tandempath = tandem_path
        self.refrowcount = 0
        self.tandemrowcount = 0
    
    def __enter__(self):
        self.reffile = open(self.refpath, 'rb')
        self.tandemfile = open(self.tandempath, 'rb')
        
        self.refreader = csv.reader(self.reffile, delimiter='\t')
        self.tandemreader = csv.reader(self.tandemfile, delimiter='\t')
        
        self.refheaders = self.refreader.next()
        self.tandemheaders = self.tandemreader.next()
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        self.reffile.close()
        self.tandemfile.close()
        
    def refnext(self):
        self.refrowcount += 1
        refrow = self.refreader.next()
        return refrow
    
    def tandemnext(self):
        self.tandemrowcount += 1
        tandemrow = self.tandemreader.next()
        return tandemrow
    
    
        
        
        
        
        