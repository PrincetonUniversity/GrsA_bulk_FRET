import os,shutil,sys

def getPDB(psf,frnum,lbnum,simdcd):
        #pull a pdb file out of a dcd simulation trajectory     
        convertCommand = "catdcd -o %d.pdb -otype pdb \
                                 -s %s -stype psf \
                                 -first %d -last %d %s > outprodcd.log" \
                                 %(lbnum,psf,frnum,frnum,simdcd)
        os.system(convertCommand)

if __name__ == "__main__":
    if os.path.exists('pdbfiles/') == False:
            os.mkdir('pdbfiles/')
    os.chdir('pdbfiles')
    #bad hard coding
    for val in range(1000):
        getPDB('../linker.psf',int(val+1),int(val+1),'../linker_ld.dcd')
