import sys
import re

""" a lame attempt for a Python representaiton of a ped file """
""" Ped file description http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped """
class PedFile(object):
    def __init__(self):
        """ family id """
        self.fid=''
        """ individual id """
        self.iid=None
        """ paternal id """
        self.pid=None
        """ maternal id """
        self.mid=None
        """ sex 1=male, 2 = female """
        self.sex=None
        """ -9 missing 0 unaffected 1 affected """
        self=pheno=None
        """ the first 6 are required, the rest are genotypes which we put in a list """
        self.genotypes=[]

        self.mother=None
        self.father=None


    def __init__(self, pedstring):
        fields = pedstring.split('\t')
        (fid, iid, pid, mid, sex, pheno) = fields[0:6]
        self.genotypes=fields[7:]
        self.fid= (fid)
        self.iid= (iid)
        self.pid= (pid)
        self.mid= (mid)
        self.sex= int(sex)
        self.pheno= int(pheno)
        
        self.mother=None
        self.father=None

    def getFid(self):
        return (self.fid)
    def getIid(self):
        return (self.iid)
    def getPid(self):
        return (self.pid)
    def getMid(self):
        return (self.mid)
    def getSex(self):
        return (self.sex)
    def getPheno(self):
        return int(self.pheno)
    def getGenotypes(self):
        return int(self.genotype)

    def setFid(self, fid):
        self.fid=fid
    def setIid(self, iid):
        self.iid=iid
    def setPid(self,pid):
        self.pid=pid
    def setMid(self,mid):
        self.mid=mid
    def setSex(self,sex):
        self.sex=sex
    def setPheno(self,pheno):
        self.pheno=pheno
    def setGenotypes(self, genotypes):
        self.genotypes=genotypes
    def setMother(self, pedobj):
        self.mother= pedobj
    def setFather(self, pedobj):
        self.father=pedobj


    def toString(self):
        outstring="\t".join( [ str(self.fid), str(self.iid), self.pid, self.mid, str(self.sex), str(self.pheno)])
        gstring="\t".join(self.genotypes)
        return outstring + "\t" + gstring

    def toStringNoPheno(self):
        outstring="\t".join( [ str(self.fid), str(self.iid), self.pid, self.mid, str(self.sex)])
        gstring="\t".join(self.genotypes)
        return outstring + "\t" + gstring


class Pedigree(object):
    def __init__(self):
        """ a Pedigree is made of founders and non-founders; non-founders have both parents in pedigree """
        self.founders=[]
        self.nonfounders=[]
        self.ped_dict={}
        self.size=0

    def __init__(self, pedfilename):
        """ given a filehandle to a ped file, initialize the Pedigree object """
        self.founders=[]
        self.nonfounders=[]
        self.ped_dict={}
        self.size=0
        fh=open(pedfilename, 'r')
        for line in fh:
            pedobject=PedFile(line.strip() )
            
            self.add(pedobject)

    def add(self, pedobj):
        """ check to see if pedobj is in the Pedigree """
        #print pedobj.getPid(), pedobj.getMid(), type(pedobj.getPid()), type(pedobj.getMid())

        if  pedobj.getMid() == '0' and pedobj.getPid() == '0' :
            sys.stderr.write("founder ...\n")
            self.founders.append(pedobj)
            self.ped_dict[ pedobj.getIid() ] = pedobj
        else:
            sys.stderr.write("nonfounder ...\n")
            self.nonfounders.append(pedobj)
            #pedobj.setMother( self.ped_dict[ pedobj.getMid()])
           # pedobj.setMother( self.ped_dict[ pedobj.getPid()])
            self.ped_dict[ pedobj.getIid() ] = pedobj


    def getiIds(self):
        founderids=map(lambda x: x.getIid(), self.founders)
        nonfounders=map(lambda x: x.getIid(), self.nonfounders)
        ids=founderids + nonfounders
        return founderids + nonfounders
    """ given and individ id, check if itsin the pedigree and set its list of genotypes """
    def addGenotypes(self, genotypes, iid):
       
        ids=self.ped_dict.keys()
        if iid not in ids:
            sys.stderr.write(iid + " not in list of iids of Pedigree object!\n")
            return
        else:
            self.ped_dict[iid].setGenotypes(genotypes)
            
    
    def dump(self):
        for s in self.founders:
            print s.toString()
        for s in self.nonfounders:
            print s.toString()

    def dumpToFile(self,fh):
        for s in self.founders:
            outstr=s.toString()
            fh.write(outstr+"\n")
        for s in self.nonfounders:
            outstr=s.toString()
            fh.write(outstr+"\n")

    def dumpToFileNoPheno(self,fh):
        for s in self.founders:
            outstr=s.toStringNoPheno()
            fh.write(outstr+"\n")
        for s in self.nonfounders:
            outstr=s.toStringNoPheno()
            fh.write(outstr+"\n")
   
        

