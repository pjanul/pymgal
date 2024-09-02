import numpy as np
from struct import unpack
from os import fstat, path
from glob import glob
import h5py
nmets = 11


def readsnapgd(filename, block, endian=None, quiet=False, longid=False, nmet=11,
                fullmass=False, mu=None, fmt=None, ptype=None, rawdata=False):
    """
    readsnapgd(filename,block,endian=None,quiet=None,longid=None,met=None, fmt=None)
        read snapshot files and new subfind files, return any block result you need.

    Parameters:
    ---------------
        filename: path plus full file name. e.g.  /your/dir/snap_009.0
        block: The block you want to read, e.g. "HEAD". Look for more info with block == "INFO"
        little endian: ">", big endian : "<", other/default : "=" or "@"
        longid: Is the particle ID saved in long long (uint64)? Default : False
        nmet: Specify how many different matels are produced in the simulation, default: 11
        fullmass: return all mass of particles inorder of saved particle position
                  False(default): return only mass block
        mu: mean_molecular_weight. Specify this value for gas temperature.
                  It will be ignored when you have NE block in your simulatin data.
        fmt: default or 1: G3 format with blocks; 0: G2 format; -1: new subfind results.
        ptype: read only specified particle type: 0: gas, 1: DM, 2: , 3: , 4: star, 5: bh
        rawdata: default False. If True, retrun the binary data in str, which need unpack yourself.

    Notes:
    ------------
    The old parameter met "z", is deprecated. If you need metal in z instead of elements,
    simply put 'Z   ' for the block.

    For these snapshots which are more than 4 Gb, i.e. the data size (bytes) indicator,
    which is just ahead of the data block, is negative, you can use `ptype=1` to overcome
    the error in reading the data.
    """

    if endian is None:
        npf = open(filename, 'rb')
        (rhead,) = unpack('<I', npf.read(4))
        npf.close()
        if (rhead == 134217728) or (rhead == 65536):
            endian = '>'
        elif (rhead == 8) or (rhead == 256):
            endian = '<'
        else:
            raise ValueError("Don't know the file endian with this value %d." % rhead)
        if not quiet:
            print("detected file format endian = ", endian)

    if fmt is None:
        # try to get the format
        npf = open(filename, 'rb')
        bs1 = unpack(endian + 'i', npf.read(4))[0]
        if bs1 == 256:
            fmt = 0
            if not quiet:
                print("Snapshot with Gadget 2 format")
        elif bs1 == 8:
            fmt = 1
            if not quiet:
                print("Snapshot with Gadget 3 format with blocks")
        else:
            print("Not knowing what is this value ", bs1, "still assuming format with block")
            print("This may have incorrect results, better have a check on endian.")
            fmt = 1
        npf.close()

    global nmets
    if nmets != nmet:
        nmets = nmet
        print('Nmetals : ', nmets)

    # read header
    npf = open(filename, 'rb')
    if fmt != 0:
        bname, bsize = read_bhead(npf)

    class rhead:
        def __init__(self, npf):
            bs1 = npf.read(4)  # size of header
            del(bs1)
            self.npart = np.zeros(6, dtype='int32')
            self.npart[:] = unpack(endian + 'i i i i i i', npf.read(4 * 6))
            self.masstbl = np.zeros(6, dtype='float64')
            self.masstbl[:] = unpack(endian + 'd d d d d d', npf.read(8 * 6))
            self.Time, self.Redshift = unpack(endian + 'd d', npf.read(2 * 8))
            self.F_Sfr, self.F_Feedback = unpack(endian + 'i i', npf.read(2 * 4))
            self.totnum = np.zeros(6, dtype='int64')
            self.totnum[:] = unpack(endian + 'i i i i i i', npf.read(6 * 4))
            self.F_Cooling, self.Numfiles = unpack(endian + 'i i', npf.read(2 * 4))
            self.Boxsize, self.Omega0, self.OmegaLambda, self.HubbleParam = unpack(endian + 'd d d d', npf.read(4 * 8))
            self.F_StellarAge, self.F_Metals = unpack(endian + 'i i', npf.read(2 * 4))
            self.nallHW = np.zeros(6, dtype='int32')
            self.nallHW[:] = unpack(endian + 'i i i i i i', npf.read(6 * 4))
            for i in range(6):
                if self.nallHW[i]>0:
                    self.totnum[i]+=self.nallHW[i]<<32
            self.F_entr_ics = unpack(endian + 'i', npf.read(4))[0]

    hd = rhead(npf)
    npf.close()

    if block == 'HEAD':
        return hd

    if block == 'IDTP':  # Particle type
        idtype = np.zeros(hd.npart.sum(), dtype=np.int32)
        nn = 0
        for i, j in enumerate(hd.npart):
            if j > 0:
                idtype[nn:nn + j] = i
                nn += j
        return(idtype)
    else:
        if fmt >= 0:
            if ptype is not None:
                if ptype == 0:
                    pty = [0, hd.npart[0]]
                elif ptype == 1:
                    pty = [hd.npart[0], hd.npart[1]]
                elif ptype == 2:
                    pty = [np.sum(hd.npart[:2]), hd.npart[2]]
                elif ptype == 3:
                    pty = [np.sum(hd.npart[:3]), hd.npart[3]]
                elif ptype == 4:
                    pty = [np.sum(hd.npart[:4]), hd.npart[4]]
                elif ptype == 5:
                    pty = [np.sum(hd.npart[:5]), hd.npart[5]]
                else:
                    raise ValueError("Don't accept ptype value %d" % ptype)
            else:
                pty = None  # the same as ptype

            if block == "MASS":
                idg0 = (hd.npart > 0) & (hd.masstbl <= 0)
                if fullmass:
                    if len(hd.npart[idg0]) == 0:  # No Mass block!
                        idg1 = (hd.npart > 0) & (hd.masstbl > 0)
                        if len(hd.npart[idg1]) == 1:
                            return hd.masstbl[idg1]
                        else:  # multi masstble
                            totmass = np.zeros(np.sum(hd.npart, dtype='int64'), dtype='float32')
                            countnm = 0
                            for i in np.arange(6):
                                if hd.npart[i] > 0:
                                    totmass[countnm:countnm + hd.npart[i]] = hd.masstbl[i]
                                    countnm += hd.npart[i]
                            return totmass
                elif ptype is not None:
                    if (hd.npart[ptype] > 0) & (hd.masstbl[ptype] > 0):
                        return hd.masstbl[ptype]
                else:
                    if len(hd.npart[idg0]) == 0:  # No Mass block!
                        return hd.masstbl

        npf = open(filename, 'rb')
        subdata = read_block(npf, block, endian, quiet, longid, fmt, pty, rawdata)
        if subdata is not None:  # we have subdata
            if block == "MASS":  # We fill the mass with the mass tbl value if needed
                npf.close()
                idg0 = (hd.npart > 0) & (hd.masstbl > 0)
                if (len(hd.npart[idg0]) > 0) and (fullmass):
                    totmass = np.zeros(np.sum(hd.npart, dtype='int64'), dtype='float32')
                    bgc = 0
                    subc = 0
                    for k in np.arange(6):
                        if hd.npart[k] > 0:
                            if(hd.masstbl[k] > 0):
                                totmass[bgc:bgc + hd.npart[k]
                                        ] = np.zeros(hd.npart[k], dtype='float32') + hd.masstbl[k]
                            else:
                                totmass[bgc:bgc + hd.npart[k]] = subdata[subc:subc + hd.npart[k]]
                                subc += hd.npart[k]
                            bgc += hd.npart[k]
                    return totmass
                else:
                    if ptype is not None:
                        if (hd.npart[ptype] == 0) or (hd.masstbl[ptype] > 0):
                            print("This is can not be! hd.npart[ptype] is ",
                                  hd.npart[ptype], "masstbl[ptype] is ", hd.masstbl[ptype])
                            print("I return 0")
                            return(None)
                        else:
                            startc = 0
                            endc = 0
                            for ii in range(ptype + 1):
                                if (hd.npart[ii] > 0) and (hd.masstbl[ii] <= 0):
                                    startc = endc
                                    endc += hd.npart[ii]
                            return(subdata[startc:endc])
                    return subdata
            elif ((block == "Z   ") or (block == "ZTOT") or (block == "Zs  ")) and (ptype is not None):
                if ptype == 0:
                    return subdata[:hd.npart[0]]
                elif ptype == 4:
                    return subdata[hd.npart[0]:]
                else:
                    raise ValueError(
                        "The given ptype %d is not accepted for metallicity block %s.", ptype, block)
            else:
                npf.close()
                return subdata
        else:  # No subdata returned
            if block == 'TEMP':  # No temperature block. Try to calculate the temperature from U
                temp = read_block(npf, "U   ", endian, 1, longid, fmt, pty, rawdata)
                if temp is None:
                    print("Can't read gas Temperature (\"TEMP\") and internal energy (\"U   \")!!")
                else:
                    xH = 0.76  # hydrogen mass-fraction
                    yhelium = (1. - xH) / (4 * xH)
                    NE = read_block(npf, "NE  ", endian, 1, longid, fmt, pty, rawdata)
                    if NE is None:
                        # we assume it is NR run with full ionized gas n_e/nH = 1 + 2*nHe/nH
                        if mu is None:
                            mean_mol_weight = (1. + 4. * yhelium) / (1. + 3 * yhelium + 1)
                        else:
                            mean_mol_weight = mu
                    else:
                        mean_mol_weight = (1. + 4. * yhelium) / (1. + yhelium + NE)
                    v_unit = 1.0e5 # * np.sqrt(hd.Time)       # (e.g. 1.0 km/sec)
                    prtn = 1.67373522381e-24  # (proton mass in g)
                    bk = 1.3806488e-16        # (Boltzman constant in CGS)
                    npf.close()
                    return(temp * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk)
            elif ((block == "Z   ") or (block == "ZTOT")):
                # no "Z   " in the data, which needs to calculate it from "Zs  " block
                subdata = read_block(npf, "Zs  ", endian, True, longid, fmt, pty, rawdata)
                if subdata is None:
                    raise ValueError("Can't find the 'Zs  ' block for calculate metallicity!")
                if ptype == 0:
                    if hd.masstbl[0] > 0:
                        mass = np.zeros(hd.npart[0], dtype=hd.masstbl.dtype) + hd.masstbl[0]
                    else:
                        mass = read_block(npf, "MASS", endian, True, longid,
                                          fmt, [0, 0], rawdata)[0:hd.npart[0]]
                    npf.close()
                    # return
                    # np.sum(subdata[0:hd.npart[0],1:],axis=1)/(mass[0:hd.npart[0]]-np.sum(subdata[0:hd.npart[0],:],axis=1))
                    # old version with z = M_z/M_H why?
                    # MASS block do not accept pty, all mass are returned!
                    return np.sum(subdata[0:hd.npart[0], 1:], axis=1) / mass
                elif ptype == 4:
                    # have to use initial mass because the metal block include SN metals.
                    im = read_block(npf, "iM  ", endian, True, longid, fmt, pty, rawdata)
                    npf.close()
                    # return
                    # np.sum(subdata[hd.npart[0]:,1:],axis=1)/(im-np.sum(subdata[hd.npart[0]:,:],axis=1))
                    # old version with z = M_z/M_H why?
                    return np.sum(subdata[hd.npart[0]:, 1:], axis=1) / im
                else:
                    zs = np.zeros(hd.npart[0] + hd.npart[4], dtype=subdata.dtype)
                    if hd.masstbl[0] > 0:
                        mass = np.zeros(hd.npart[0], dtype=hd.masstbl.dtype) + hd.masstbl[0]
                    else:
                        mass = read_block(npf, "MASS", endian, True, longid,
                                          fmt, [0, 0], rawdata)[0:hd.npart[0]]
                    # zs[0:hd.npart[0]]=np.sum(subdata[0:hd.npart[0],1:],axis=1)/(mass[0:hd.npart[0]]-np.sum(subdata[0:hd.npart[0],:],axis=1))
                    # old version
                    zs[0:hd.npart[0]] = np.sum(subdata[0:hd.npart[0], 1:], axis=1) / mass

                    im = read_block(npf, "iM  ", endian, True, longid, fmt, pty, rawdata)
                    # zs[hd.npart[0]:]=np.sum(subdata[hd.npart[0]:,1:],axis=1)/(im-np.sum(subdata[hd.npart[0]:,:],axis=1))
                    zs[hd.npart[0]:] = np.sum(subdata[hd.npart[0]:, 1:], axis=1) / im
                    mass, im, subdata = 0, 0, 0
                    npf.close()
                    return zs

            if not quiet:
                print("No such blocks!!! or Not add in this reading!!!", block)
            npf.close()
            return(None)


# Read Block
def read_block(npf, block, endian, quiet, longid, fmt, pty, rawdata):
    from readsnapsgl import readsnapgd # import rhead
    global nmets
    endf = fstat(npf.fileno()).st_size

    header=readsnapgd(npf.name,'HEAD')
    bname = 'BLOCK_NAME'
    if fmt == 0:
        npf.seek(8 + 256)  # skip block(16) + header (264)
    elif fmt == 1:
        npf.seek(16 + 8 + 256)  # skip header (264)
    loopnum = 0
    # while bname!='EOFL' :   #Ending block
    while npf.tell() < endf:  # End of file
        if fmt != 0:
            bname, bsize = read_bhead(npf)
            bsize = npf.read(4)
            bsize = unpack(endian + 'i', bsize)[0]
            npf.seek(npf.tell() - 4)
        else:
            bsize = npf.read(4)
            bsize = unpack(endian + 'i', bsize)[0]
            npf.seek(npf.tell() - 4)

            if (block == 'POS ') and (loopnum == 0):
                return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
            elif (block == 'VEL ') and (loopnum == 1):
                return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
            elif (block == 'ID  ') and (loopnum == 2):
                if longid:
                    return read_bdata(npf, 1, np.dtype('uint64'), endian, pty)
                else:
                    return read_bdata(npf, 1, np.dtype('uint32'), endian, pty)
            elif (block == 'MASS') and (loopnum == 3):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'U   ') and (loopnum == 4):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'RHO ') and (loopnum == 5):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'NE  ') and (loopnum == 6):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'NH  ') and (loopnum == 7):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'HSML') and (loopnum == 8):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'SFR ') and (loopnum == 9):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'DT  ') and (loopnum == 10):  # delayed time
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif (block == 'AGE ') and (loopnum == 11):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif ((block == 'Z   ') or (block == 'Zs  ')) and (loopnum == 12):
                if nmets <= 0:
                    return None
                elif nmets == 1:  # suppose to be metallicity z
                    return read_bdata(npf, 1, np.dtype('float32'), endian, pty)
                elif block == 'Zs  ':
                    return read_bdata(npf, nmets, np.dtype('float32'), endian, pty)
                else:
                    zs = read_bdata(npf, nmets, np.dtype('float32'), endian, pty)
                    npf.seek(4)
                    npart = unpack(endian + 'i i i i i i', npf.read(4 * 6))
                    npf.seek(264)
                    for i in range(3):
                        bs1 = unpack(endian + 'i', npf.read(4))[0]
                        npf.seek(npf.tell()+bs1+4)
                    # bs1 = unpack(endian + 'i', npf.read(4))[0]
                    mass = read_bdata(npf, 1, np.dtype('float32'), endian, pty)
                    # note this only return the gas metallicity!!!
                    return np.sum(zs, axis=1)/mass
                    # return np.sum(zs[:npart[0]], axis=1)/mass[:npart[0]]
            elif loopnum > 12:
                return None
            loopnum += 1

        if bsize<=0 or header.npart.sum()*4*3>2147483647:
            print("!!!# WARNING: The read skip may not work correct here!!! bsize=", bsize, "total part number= header.npart.sum()")
            print("!!!# WARNING: Check your results!!! I am trying to use particle number with specified ptype!")
            if loopnum<3:  #Note the loopnum +=1 before, pos: loopnum=1
                bsize=12*header.npart.sum()
            elif loopnum==3:  # IDs: loopnum=3
                if longid:
                    bsize=8*header.npart.sum()
                else:
                    bsize=4*header.npart.sum()
            else:
                bsize=4*header.npart.sum()
        if not quiet:
            if fmt != 0:
                print(bname, bsize)
            else:
                print("Format 0, reading block ", block,  "skiping", bsize)

        # For reading snapshot files###
        if rawdata:
            if bname == block:
                return npf.read(unpack(endian + 'i', npf.read(4))[0])

        if bname == block == 'POS ':
            return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
        elif bname == block == 'VEL ':
            return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
        elif bname == block == 'ID  ':
            if longid:
                return read_bdata(npf, 1, np.dtype('uint64'), endian, pty)
            else:
                return read_bdata(npf, 1, np.dtype('uint32'), endian, pty)
        elif bname == block == 'MASS':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'RHO ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'NE  ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'NH  ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'SFR ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'AGE ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'POT ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'iM  ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'Zs  ':
            return read_bdata(npf, nmets, np.dtype('float32'), endian)
        elif bname == block == 'HOTT':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'Z   ':  # specified block, which saves the metallicity z
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'ZTOT':  # specified block, which saves the metallicity z
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'CLDX':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'MHI ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'TEMP':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'HSML':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'PTYP':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        # Internal energy###
        elif bname == block == 'U   ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        # INFO ####
        elif bname == block == 'INFO':  # This is print out not return in array
            bs1 = unpack(endian + 'i', npf.read(4))[0]
            buf = npf.read(bs1)
            print("Block   DataType   dim  Type0 Type1 Type2 Type3 Type4 Type5")
            cc = 0
            while cc < bs1:
                print(unpack(endian + '4s 8s i i i i i i i', buf[cc:cc + 40]))
                cc += 40
            return(1)

        # For reading new subfind files###
        elif bname == block == 'GLEN':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'GOFF':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'MTOT':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'GPOS':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'MVIR':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'RVIR':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'M25K':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'R25K':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'M500':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'R500':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'MGAS':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'MSTR':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'TGAS':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'LGAS':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'NCON':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'MCON':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'BGPO':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'BGMA':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'BGRA':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'NSUB':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'FSUB':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'SLEN':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'SOFF':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'SSUB':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'MSUB':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'SPOS':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'SVEL':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'SCM ':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'SPIN':
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif bname == block == 'DSUB':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'VMAX':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'RMAX':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'MBID':
            if longid:
                return read_bdata(npf, 1, np.dtype('int64'), endian)
            else:
                return read_bdata(npf, 1, np.dtype('uint32'), endian)
        elif bname == block == 'GRNR':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif bname == block == 'SMST':
            return read_bdata(npf, 6, np.dtype('float32'), endian)
        elif bname == block == 'SLUM':
            return read_bdata(npf, 12, np.dtype('float32'), endian)
        elif bname == block == 'SLAT':
            return read_bdata(npf, 12, np.dtype('float32'), endian)
        elif bname == block == 'SLOB':
            return read_bdata(npf, 12, np.dtype('float32'), endian)
        elif bname == block == 'DUST':
            return read_bdata(npf, 11, np.dtype('float32'), endian)
        elif bname == block == 'SAGE':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'SZ  ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'SSFR':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'PID ':
            if longid:
                return read_bdata(npf, 1, np.dtype('int64'), endian)
            else:
                return read_bdata(npf, 1, np.dtype('uint32'), endian)

        # For reading my PIAO outputs ###
        elif (bname == block == 'GOFF') or (bname == block == 'GHED') or \
            (bname == block == 'GSBL') or (bname == block == 'GSBO') or \
                (bname == block == 'SBLN') or (bname == block == 'SBOF'):
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        elif (bname == block == 'GMAS') or (bname == block == 'GRAD') or (bname == block == 'SBMS'):
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif (bname == block == 'GDCP') or (bname == block == 'GMCP') or (bname == block == 'SBPS'):
            return read_bdata(npf, 3, np.dtype('float32'), endian)
        elif (bname == block == 'GCID') or (bname == block == 'GIDS') or (bname == block == 'SBCI'):
            if longid:
                return read_bdata(npf, 1, np.dtype('int64'), endian)
            else:
                return read_bdata(npf, 1, np.dtype('uint32'), endian)

        else:
            # if fmt != 0:
            #    bsize=unpack(endian+'i',npf.read(4))[0]
            #    npf.seek(bsize+npf.tell()+4)
            # else:
            npf.seek(bsize + 8 + npf.tell())

    return None


# Read Block Head
def read_bhead(npf):
    dummy = npf.read(4)  # dummy
    if dummy == '':
        bname = "EOFL"
        bsize = '0'
    else:
        bname = npf.read(4).decode('ascii')  # label
        bsize = npf.read(4)  # size
        if npf.read(4) != dummy:
            print("header part not consistent!!")
    return bname, bsize


# Read Block data
def read_bdata(npf, column, dt, endian, pty=None):
    bs1 = unpack(endian + 'i', npf.read(4))[0]

    if pty is None:
        buf = npf.read(bs1)
    else:
        npf.seek(npf.tell() + pty[0] * dt.itemsize * column)
        buf = npf.read(pty[1] * dt.itemsize * column)
        bs1 = pty[1] * dt.itemsize * column

    if column == 1:
        arr = np.ndarray(shape=np.int32(bs1 / dt.itemsize), dtype=dt, buffer=buf)
    else:
        arr = np.ndarray(shape=(np.int32(bs1 / dt.itemsize / column), column),
                         dtype=dt, buffer=buf)

    if (endian == '=') or (endian == '<'):  # = and < gives the same result
        return arr
    else:
        return arr.byteswap()


def readhdf5head(filename, quiet=False):
    if not quiet:
        print('Reading %s file with Header' % filename)
    fo = h5py.File(filename, 'r')
    class rhead:
        def __init__(self, npf):
            self.npart = npf['Header'].attrs['NumPart_ThisFile']
            self.masstbl = npf['Header'].attrs['MassTable']
            self.Time = npf['Header'].attrs['Time']
            self.Redshift = npf['Header'].attrs['Redshift']
            self.F_Sfr = npf['Header'].attrs['Flag_Sfr']
            self.F_Feedback = npf['Header'].attrs['Flag_Feedback']
            self.totnum = npf['Header'].attrs['NumPart_Total']
            self.F_Cooling = npf['Header'].attrs['Flag_Cooling']
            self.Numfiles = npf['Header'].attrs['NumFilesPerSnapshot']
            self.Boxsize = npf['Header'].attrs['BoxSize']
            self.Omega0 = npf['Header'].attrs['Omega0']
            self.OmegaLambda = npf['Header'].attrs['OmegaLambda']
            self.HubbleParam = npf['Header'].attrs['HubbleParam']
            self.F_StellarAge = npf['Header'].attrs['Flag_StellarAge']
            self.F_Metals = npf['Header'].attrs['Flag_Metals']
            self.F_DoublePrecision = npf['Header'].attrs['Flag_DoublePrecision']
    hd = rhead(fo)
    fo.close()
    return hd

def readhdf5data(filename, block, quiet=False, ptype=None):
    if not quiet:
        print('Reading file ', filename, ' with data block ', block,' for type ', ptype)
    fo = h5py.File(filename, 'r')

    if 'PartType0' in fo.keys():
        if (block.lower() == 'temperature') & ('temperature' not in [k.lower() for k in fo['PartType0'].keys()]):
            if 'InternalEnergy' not in fo['PartType0'].keys():
                print("Can't read gas Temperature and internal energy as both are not in gas properties: !!", fo['PartType0'].keys())
                fo.close()
                return None
            else:
                temp = fo['PartType0/InternalEnergy'][:]
                if 'Metallicity' not in fo['PartType0'].keys():
                    xH = 0.76  # hydrogen mass-fraction
                else:
                    Metals=fo['PartType0/Metallicity'][:]
                    if Metals.shape[1] == 11:
                        xH = 1-Metals[:,0]-Metals[:,1]
                    elif Metals.shape[1] == 34:
                        if not quiet: print('Using Chem5 model in SIMBA!')
                        xH = Metals[:,1]
                yhelium = (1. - xH) / (4 * xH)
                if 'ElectronAbundance' not in fo['PartType0'].keys(): # we assume it is NR run with full ionized gas n_e/nH = 1 + 2*nHe/nH
                    mean_mol_weight = (1. + 4. * yhelium) / (1. + 3 * yhelium + 1)
                else:
                    NE = fo['PartType0/ElectronAbundance'][:]
                    mean_mol_weight = (1. + 4. * yhelium) / (1. + yhelium + NE)
                v_unit = 1.0e5 # * np.sqrt(fo['Header'].attrs['Time'])       # (e.g. 1.0 km/sec)
                prtn = 1.67373522381e-24  # (proton mass in g)
                bk = 1.3806488e-16        # (Boltzman constant in CGS)
                fo.close()
                return(temp * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk)

    if isinstance(ptype, type(0)):
        if 'PartType'+str(ptype) in fo.keys():
            if block in fo['PartType'+str(ptype)].keys():
                data = fo['PartType'+str(ptype)+'/'+block][:]
            else:
                print(block, ' is not in PartType'+str(ptype), '!!')
                fo.close()
                return None
        else:
            print('PartType'+str(ptype), 'is not in this HDF5 file: ', filename, ' !!')
            fo.close()
            return None
    else:
        if ptype is None:  #read all types of data in the file
            PT=list(fo.keys())
            PT.remove('Header')
        else:
            PT=[]
            for i in ptype:
                PT.append('PartType'+str(ptype))

        for i, ptn in enumerate(PT):
            if i == 0:
                if ptn in fo.keys():
                    if block in fo[ptn].keys():
                        data = fo[ptn+'/'+block][:]
                    else:
                        print('# WARNING: ', block, ' is not in ', ptn, '!!')
                else:
                    print('# WARNING: ', ptn, ' is not in this HDF5 file: ', filename, '!!')
            else:
                if ptn in fo.keys():
                    if block in fo[ptn].keys():
                        data = np.append(data, fo[ptn+'/'+block][:], axis=0)
                    else:
                        print('# WARNING: ', block, ' is not in ', ptn, '!!')
                else:
                    print('# WARNING: ', ptn, ' is not in this HDF5 file: ', filename, '!!')

    fo.close()
    return data

def cut_sim_region(spos, center, radius, boxsize, restrict_r):
    for i in range(3):
        if center[i]+radius > boxsize and center[i] - radius < 0:
            raise ValueError('Too large region!! center[i]+radius: %f, center[i] - radius: %f, Please use a smaller raidus' % (center[i]+radius, center[i] - radius))
        if center[i]+radius > boxsize:
            ids=spos[:,i] <= center[i]+radius - boxsize
            spos[ids,i] += boxsize
        if center[i] - radius < 0:
            ids = spos[:,i] >= center[i] - radius + boxsize
            spos[ids,i] -= boxsize
    if restrict_r:
        return np.sqrt(np.sum((spos - center)**2, axis=1)) <= radius
    else:
        # ids = r <= self.radius*np.sqrt(2)  # increase to get all the projected data
        # Now using cubic box to get the data
        return  (spos[:, 0] > center[0] - radius) & \
                (spos[:, 0] <= center[0] + radius) & \
                (spos[:, 1] > center[1] - radius) & \
                (spos[:, 1] <= center[1] + radius) & \
                (spos[:, 2] > center[2] - radius) & \
                (spos[:, 2] <= center[2] + radius)


# read all snapshots
def readsnap(filename, block, endian=None, quiet=True, longid=False, nmet=11, totnum=None,
             fullmass=False, mu=None, fmt=None, ptype=None, rawdata=False, randf=None,
             center=None, radius=None, restrict_r=True):
    """
    readsnap(filename, block, endian=None, quiet=False, longid=False, nmet=11,
             fullmass=False, mu=None, fmt=None, ptype=None, rawdata=False
        read multiple snapshot files and new subfind files, return any block data in whole simulation.

    Parameters:
    ---------------
        filename: path plus full file name. e.g.  /your/dir/snap_009.0
        block: The block you want to read, e.g. "HEAD". Look for more info with block == "INFO"
        little endian: ">", big endian : "<", other/default : "=" or "@"
        longid: Is the particle ID saved in long long (uint64)? Default : False
        nmet: Specify how many different matels are produced in the simulation, default: 11
        totnum: If specified, this number instead of the one coming from header (even with ptype!=None)
                is used to create the array.
        fullmass: return all mass of particles inorder of saved particle position
                  False(default): return only mass block
        mu: mean_molecular_weight. Specify this value for gas temperature.
                  It will be ignored when you have NE block in your simulatin data.
        fmt: default or 1: G3 format with blocks; 0: G2 format; -1: new subfind results.
        ptype: read only specified particle type: 0: gas, 1: DM, 2: , 3: , 4: star, 5: bh
        rawdata: default False. If True, retrun the binary data in str, which need unpack yourself.
        randf: float (0-1], a random fraction of the total data will be returned.
        -------------------------------------------------------------------
        Note these options are not fully compile with randf. No postion data can also be queried.
        center      : The center of a sphere for the data you want to get.
                      Default : None, whole data will load.
        radius      : The radius of a sphere for the data you want to get.
                      Default : None, whole data will be used.
        restrict_r  : Using the exact radius to cut out simulation data.
                      Default : True. Otherwise, a cubic (2 *radius) data is cut out to fill the output fits image.

    Notes:
    ------------
    The old parameter met "z", is deprecated. If you need metal in z instead of elements,
    simply put 'Z   ' for the block.

    For these snapshots which are more than 4 Gb, i.e. the data size (bytes) indicator,
    which is just ahead of the data block, is negative, you can use `ptype=1` to overcome
    the error in reading the data.

    For the gadget2 snapshots files, please check the order of the reading is correct or not (line 296-324) for you data.
    """

    if path.isfile(filename): ## only one simulation file
        filenum=1
        filename=[filename]
    else:
        filename=glob(filename+'*')
        filenum=len(filename)
        if len(filename) > 1:
            if '.hdf5' in ','.join(filename).lower():
                filename = [ x for x in filename if x[-4:].lower() == 'hdf5']  #exclude the other files
                filenum=len(filename)
            elif '.0' in ','.join(filename).lower():
                filename = [ x for x in filename if x.split('.')[-1].isnumeric()]  #exclude the other files
                filenum=len(filename)
        else:
            raise ValueError("Can not find file:", filename,filename+"*")

    if not quiet:
        print('reading files: ', filename)

    if filename[0][-4:].lower() == 'hdf5':
        head = readhdf5head(filename[0], quiet=quiet)
        if block == 'HEAD':
            return head
        elif (block == 'IDTP') or (block == 'IDTypes'):
            data=np.zeros(head.npart[0],dtype=np.int32)
            for i in np.arange(1,6):
                if head.npart[i]>0:
                    data = np.append(data, np.ones(head.npart[i],dtype=np.int32)*i, axis=0)
            return data

    else:
        head=readsnapgd(filename[0], 'HEAD', endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                         fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
        if block == 'HEAD':
            return head

    if filenum !=1 and head.Numfiles != filenum: # allow to read only one file
        raise ValueError("The number of files (%i) do not fit to the one in snapshot header (%i) !! please check!" % (filenum, head.Numfiles))

    if totnum is None:
        if ptype is None:
            totnum = head.totnum.sum()
        else:
            totnum = head.totnum[ptype]

    if randf is not None:
        Rtotnum = np.int64(totnum * randf)
    else:
        Rtotnum = totnum

    cut_region=False
    if (center is not None) and (radius is not None):
        cut_region=True

    if filenum == 1: # only one file
        if randf is not None:
            rids=np.random.random_integers(0, high=totnum-1, size=Rtotnum)

        if filename[0][-4:].lower() == 'hdf5':
            if randf is None:
                if cut_region:
                    pos = readhdf5data(filename[0], 'Coordinates', quiet=quiet, ptype=ptype)
                    cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                    if block == 'Coordinates':
                        return pos[cids]
                    else:
                        return readhdf5data(filename[0], block, quiet=quiet, ptype=ptype)[cids]
                else:
                    return readhdf5data(filename[0], block, quiet=quiet, ptype=ptype)
            else:
                return readhdf5data(filename[0], block, quiet=quiet, ptype=ptype)[rids]
        else:
            if randf is None:
                if cut_region:
                    pos = readsnapgd(filename[0], 'POS ', endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                       fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
                    cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                    if block == 'POS ':
                        return pos[cids]
                    else:
                        return readsnapgd(filename[0], block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                           fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)[cids]
                else:
                    return readsnapgd(filename[0], block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                       fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
            else:
                return readsnapgd(filename[0], block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                   fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)[rids]

    else:  # multiple snapshot names
        if ptype is not None:
            if head.totnum[ptype]==0:
                print('snapshot does not containt this type %d particles' % ptype)
                return None

        countN,countF=0,0
        tmp=None
        if filename[0][-4:].lower() == 'hdf5':
            while tmp is None:
                if cut_region:
                    pos = readhdf5data(filename[countF], 'Coordinates', quiet=quiet, ptype=ptype)
                    cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                    if block == 'Coordinates':
                        tmp = pos[cids]
                    else:
                        tmp = readhdf5data(filename[countF], block, quiet=quiet, ptype=ptype)[cids]
                else:
                    tmp = readhdf5data(filename[countF], block, quiet=quiet, ptype=ptype)
                countF+=1
        else:
            while tmp is None:
                if cut_region:
                    pos = readsnapgd(filename[countF], 'POS ', endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                       fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
                    cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                    if block == 'POS ':
                        tmp = pos[cids]
                    else:
                        tmp = readsnapgd(filename[countF], block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                   fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)[cids]
                else:
                    tmp = readsnapgd(filename[countF], block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                   fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
                countF+=1

        if len(tmp.shape)==1:
            if randf is not None:
                data=np.zeros(Rtotnum,dtype=tmp.dtype)
            else:
                if cut_region: # don't need all data
                    if radius*4<head.Boxsize:
                        data=np.zeros(np.int64(head.totnum.sum()*(radius*4/head.Boxsize)**3),dtype=tmp.dtype)
                    else:
                        data=np.zeros(totnum,dtype=tmp.dtype)
                else:
                    data=np.zeros(totnum,dtype=tmp.dtype)
        elif len(tmp.shape)==2:
            if randf is not None:
                data=np.zeros((Rtotnum,tmp.shape[1]),dtype=tmp.dtype)
            else:
                if cut_region: # don't need all data
                    if radius*4<head.Boxsize:
                        data=np.zeros((np.int64(head.totnum.sum()*(radius*4/head.Boxsize)**3), tmp.shape[1]),dtype=tmp.dtype)
                    else:
                        data=np.zeros((totnum, tmp.shape[1]),dtype=tmp.dtype)
                else:
                    data=np.zeros((totnum, tmp.shape[1]),dtype=tmp.dtype)

        if randf is not None:
            rN = np.int64(tmp.shape[0]*randf)
            rids=np.random.random_integers(0, high=tmp.shape[0]-1, size=rN)
            data[countN:countN+rN] = tmp[rids]
            countN+=rN
        else:
            data[countN:countN+tmp.shape[0]] = tmp
            countN+=tmp.shape[0]

        for i,fbase in enumerate(filename):
            if i >= countF: # new reading from previous checking
                if fbase[-4:].lower() == 'hdf5':
                    if cut_region:
                        pos = readhdf5data(fbase, 'Coordinates', quiet=quiet, ptype=ptype)
                        cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                        if block == 'Coordinates':
                            tmp = pos[cids]
                        else:
                            tmp = readhdf5data(fbase, block, quiet=quiet, ptype=ptype)[cids]
                    else:
                        tmp = readhdf5data(fbase, block, quiet=quiet, ptype=ptype)
                else:
                    if cut_region:
                        pos = readsnapgd(fbase, 'POS ', endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                           fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
                        cids= cut_sim_region(pos, center, radius, head.Boxsize, restrict_r)
                        if block == 'POS ':
                            tmp = pos[cids]
                        else:
                            tmp = readsnapgd(fbase, block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                       fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)[cids]
                    else:
                        tmp = readsnapgd(fbase, block, endian=endian, quiet=quiet, longid=longid, nmet=nmet,
                                          fullmass=fullmass, mu=mu, fmt=fmt, ptype=ptype, rawdata=rawdata)
                if tmp is not None:
                    if randf is not None:
                        rN = np.int64(tmp.shape[0]*randf)+1
                        if rN+countN > Rtotnum:  # make sure not exceeding the random array size
                            rN = Rtotnum - countN
                        if rN <= 0:
                            print('No data is read for ', i, ' because random size =', rN)
                        rids=np.random.random_integers(0, high=tmp.shape[0]-1, size=rN)
                        data[countN:countN+rN] = tmp[rids]
                        countN+=rN
                    else:
                        data[countN:countN+tmp.shape[0]] = tmp
                        countN+=tmp.shape[0]

        return(data[:countN]) # make sure no 0 in the final result
