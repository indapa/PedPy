class Ped(object):
    """ represents a pedigree record in a Ped file http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped """

    def __init__(self,famid='.', indv='.', paternal='.', maternal='.', sex='.', phenotype='.'):
        """ attributes of a Ped object """
        self.famid=famid
        self.individ=indv
        self.pid=paternal
        self.mid=maternal
        self.sex=sex
        self.pheno=phenotype


    def setfamid(self,famid):
        self.famid=famid

    def setindvid(self,indv):
        self.individ=indv

    def setmid(self,mid):
         self.mid=mid

    def setpid(self,pid):
         self.pid=pid

    def setsex(self,sex):
         self.sex=sex

    def setpheno(self,pheno):
         self.pheno=pheno



    def getfamid(self):
        return self.famid

    def getid(self):
        return self.individ

    def getmid(self):
        return self.mid

    def getpid(self):
        return self.pid

    def getsex(self):
        return self.sex

    def getpheno(self):
        return self.pheno


    def isFounder(self):
        return self.pid == '0' and self.mid == '0'

    def getParents(self):

        return ( self.pid, self.mid)


    def __str__(self):
        return "\t".join( [ self.famid, self.individ, self.pid, self.mid, self.sex, self.pheno] )


class Pedfile(object):
    """ a Pedfile object has a list of Ped  objects """

    def __init__(self,filename):
        self.filename=filename
        self.fh=open(self.filename, 'r')
        self.pedlist=[]

    def parsePedfile(self):
        """ given a filehandle to a *.ped file read its contents and populate the list pedlist with Ped objects """
        for line in self.fh:
            fields=line.strip().split('\t')
            (famid, indv, pid, mid, sex, pheno)=fields[0:6]
            self.pedlist.append( Ped(famid, indv, pid, mid, sex, pheno) )

    
    def returnFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        founders=[]

        for pedobj in self.pedlist:
            if pedobj.getpid() == "0":
                founders.append(pedobj)

        return founders

    def returnFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        founderids=[]
        for pedobj in self.pedlist:
            if pedobj.getpid() == "0":
                founderids.append( pedobj.getid() )
        return founderids

    def returnNonFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        nonfounderids=[]
        for pedobj in self.pedlist:
            if pedobj.getpid() != "0":
                nonfounderids.append( pedobj.getid() )
        return nonfounderids

    def returnNonFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        nonfounders=[]

        for pedobj in self.pedlist:
            if pedobj.getpid() != "0":
                nonfounders.append(pedobj)

        return nonfounders

    def returnIndivids(self):
        """ return a list of indvi ids from the ped file """
        #samplelist=[]
        return [ ped.getid() for ped in self.pedlist]

    def getPedList(self):
        return self.pedlist
    def getTotalSize(self):
        return len(self.pedlist)

    def yieldMembers(self):
        for pedobj in self.pedlist:
            yield pedobj

    def __str__(self):
        return "\n".join( [ x.__str__() for x in self.pedlist ] )
