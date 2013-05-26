'''
Created on 2013-05-25

@author: jyeung
'''


import csv

        
class data(object):
    '''
    Reads from genome info file and tandem file.  
    '''
    
    def __init__(self, ref_path, tandem_path, tandem_output_path):
        '''
        Filepath
        '''
        self.refpath = ref_path
        self.tandempath = tandem_path
        self.outputpath = tandem_output_path
        self.refrowcount = 0
        self.tandemrowcount = 0
        self.outputrowcount = 0
    
    def __enter__(self):
        self.reffile = open(self.refpath, 'rb')
        self.tandemfile = open(self.tandempath, 'rb')
        self.outfile = open(self.outputpath, 'wb')
        
        self.refreader = csv.reader(self.reffile, delimiter='\t')
        self.tandemreader = csv.reader(self.tandemfile, delimiter='\t')
        self.outwriter = csv.writer(self.outfile, delimiter='\t')
        
        self.refheaders = self.refreader.next()
        self.tandemheaders = self.tandemreader.next()
        # Initialize output header in writecolnames() function.
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        self.reffile.close()
        self.tandemfile.close()
        self.outfile.close()
        
    def writecolnames(self, *colnames):
        header = self.tandemheaders
        for cname in colnames:
            header.append(cname)
        self.outheaders = header
        self.outwriter.writerow(header)
        
    def refnext(self):
        self.refrowcount += 1
        refrow = self.refreader.next()
        return refrow
    
    def tandemnext(self):
        self.tandemrowcount += 1
        tandemrow = self.tandemreader.next()
        return tandemrow
    
    def writerow(self, row):
        self.outputrowcount += 1
        self.outwriter.writerow(row)
    
    
        
        
        
        
        