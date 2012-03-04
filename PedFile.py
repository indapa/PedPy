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
        gstring="\t".join(fields[6::])
        # we remove any phasing/unphasing chars in the genotypes
        # . missing genotypes are denoted -1
        gstring=gstring.replace('.', '-1 -1')
        gstring=gstring.replace('/', ' ')
        gstring=gstring.replace('|', ' ')
        
        self.genotypes=gstring.split('\t')
        self.fid= (fid)
        self.iid= (iid)
        self.pid= (pid)
        self.mid= (mid)
        self.sex= int(sex)
        self.pheno= int(pheno)
        
        self.mother=None
        self.father=None



    def __init__(self, pedstring, nopheno=True):
        fields = pedstring.split('\t')
        
        
        (fid, iid, pid, mid, sex) = fields[0:5]
        gstring="\t".join(fields[5::])
        
        # we remove any phasing/unphasing chars in the genotypes
        # . missing genotypes are denoted -1
        gstring=gstring.replace('.', '-1 -1')
        gstring=gstring.replace('/', ' ')
        gstring=gstring.replace('|', ' ')

        self.genotypes=gstring.split('\t')
        self.fid= (fid)
        self.iid= (iid)
        self.pid= (pid)
        self.mid= (mid)
        self.sex= int(sex)
        self.pheno= -1

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
        return self.genotypes

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

#############################
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


    def __init__(self, pedfilename, mendelsoft=True):
        """ given a filehandle to a mendelsoft file, initialize the Pedigree object """
        self.founders=[]
        self.nonfounders=[]
        self.ped_dict={}
        self.size=0
        fh=open(pedfilename, 'r')
        for line in fh:
            nopheno=True
            pedobject=PedFile(line.strip(), nopheno=True )

            self.add(pedobject)

    """ return a list of affecteds and unaffecteds from the Pedigree """
    def getAffecteds(self):
        affecteds=[]
        unaffecteds=[]
        for key in self.ped_dict.keys():
            if self.ped_dict[key].getPheno() == 1:
                unaffecteds.append(key)
            elif self.ped_dict[key].getPheno() == 2:
                affecteds.append(key)
            else:
                pass
        return (affecteds, unaffecteds)
    
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

    """ get matrix of genotypes; i.e. a list of lists """

    def getGenotypeMatrix(self):
        genotype_matrix=[]
        for founder in self.founders:
           
            pedobj=self.ped_dict[founder.getIid()]
            curr_genotypes=pedobj.getGenotypes()
           
            curr_genotypes.insert(0,founder.getIid())
            genotype_matrix.append(curr_genotypes)
        for nonfounder in self.nonfounders:
           
            pedobj=self.ped_dict[nonfounder.getIid()]
            curr_genotypes=pedobj.getGenotypes()
            
            curr_genotypes.insert(0,nonfounder.getIid())
            genotype_matrix.append(curr_genotypes)

        return genotype_matrix

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
   
        
""" a lame attempt for a Python representaiton of a map file """
""" Map file description http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map """
class Map(object):
    def __init__(self):
        """chromosome """
        self.chrom=''
        """ marker id """
        self.id=None
        """ genetic map position in morgans"""
        self.mapos=None
        """ physical position (1-based) """
        self.pos=None

    def __init__(self, mapstring):
        fields = mapstring.split('\t')
        (chrom, id, mapos, pos) = fields[0:4]
        self.chrom=chrom
        self.id=id
        self.mapos=float(mapos)
        self.pos=int(pos)

    def getChrom(self): return self.chrom
    def getPos(self): return self.pos
    def getMap(self): return self.mapos
    def getId(self): return id
    
    def setChrom(self, chr): self.chrom=chr
    def setId(self, id): self.id=id
    def setMap(self, map): self.map=float(map)
    def setPos(self, pos): self.pos=int(pos)

    
    def toString(self):
        return "\t".join([self.chrom, self.id, str(self.mapos), str(self.pos)])

    """ chrom pos, id are represented in a map """
    def toStringVcf(self):
         
        return "\t".join([self.chrom, str(self.pos), self.id])