import numpy as np
from struct import unpack
# from time import time
from os import fstat
nmets = 11


def readsnapsgl(filename, block, endian=None, quiet=None, longid=None, met=None, nmet=11, fullmass=False, rhb=True, fmt=None, ptype=None, rawdata=False):
    """
    readsnapsgl(filename,block,endian=None,quiet=None,longid=None,met=None, fmt=None)
        read snapshot files and new subfind files, return any block result you need.
    Parameters:
        filename: path plus full file name. e.g.  /your/dir/snap_009.0
        block : The block you want to read, e.g. "HEAD". Look for more info with block == "INFO"
    #little endian : ">", big endian : "<", other/default : "=" or "@"
    #longid False or leave it empty for not using. =True if long id is needed
    #met, "z", if you need metal in z instead of elements.
    #nmet,  Specify how many different matels are produced in the simulation, default: 11
    #fullmass, True: return all mass of particles inorder of saved particle position
    #          False(default): return only mass block
    #rhb, short for return head brief. True(deflaut): only return useful head information, else: all
    #fmt, default or 1: G3 format with blocks; 0: G2 format for snap shot
    #      -1: new subfind results.
    #ptype, read only specified particle type: 0: gas, 1: DM, 2: , 3: , 4: star, 5: bh
    #rawdata, default false. If True, will retrun the read binary data in str. You need unpack yourself.
    """

    if endian is None:
        endian = "="
    if quiet is None:
        quiet = False
    else:
        quiet = True
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
            print("Not knowing what is this value ", bs1,
                  "still assuming format with block")
            print("This may have incorrect results, better have a check on endian.")
            fmt = 1
        npf.close()
    if longid is None:
        longid = False

    global nmets
    if nmets != nmet:
        nmets = nmet

    if block == 'HEAD':
        npf = open(filename, 'rb')
        if fmt != 0:
            bname, bsize = read_bhead(npf)
        bs1 = npf.read(4)  # size of header
        npart = np.zeros(6, dtype='int32')
        npart[:] = unpack(endian + 'i i i i i i', npf.read(4 * 6))
        masstbl = np.zeros(6, dtype='float64')
        masstbl[:] = unpack(endian + 'd d d d d d', npf.read(8 * 6))
        time, red = unpack(endian + 'd d', npf.read(2 * 8))
        F_sfr, F_fb = unpack(endian + 'i i', npf.read(2 * 4))
        Totnum = np.zeros(6, dtype='int64')
        Totnum[:] = unpack(endian + 'i i i i i i', npf.read(6 * 4))
        F_cool, Numfiles = unpack(endian + 'i i', npf.read(2 * 4))
        Boxsize, Omega0, OmegaLambda, Hubbleparam = unpack(
            endian + 'd d d d', npf.read(4 * 8))
        F_agn, F_metal = unpack(endian + 'i i', npf.read(2 * 4))
        NallHW = np.zeros(6, dtype='int32')
        NallHW[:] = unpack(endian + 'i i i i i i', npf.read(6 * 4))
        F_entr_ics = unpack(endian + 'i', npf.read(4))[0]
        npf.close()
        if rhb:
            return(npart, masstbl, time, red, Totnum, Boxsize, Omega0, OmegaLambda, Hubbleparam)
        else:
            return(npart, masstbl, time, red, F_sfr, F_fb, Totnum, F_cool, Numfiles, Boxsize, Omega0, OmegaLambda, Hubbleparam, F_agn, F_metal, NallHW, F_entr_ics)
    # elif (block == 'INFO') & ((fmt == -1) | (fmt == 1)):
    #     npf=open(filename,'rb')
    #     npf.seek(16+264)
    #     bname,bsize=read_bhead(npf)
    #     if bname != 'INFO':
    #         print("Not reading INFO block, read block is ", bname)
    #         print("Now, I qiut")
    #         npf.close()
    #         return(0)

    #     bs1=unpack(endian+'i',npf.read(4))[0]
    #     buf=npf.read(bs1)
    #     print("Block   DataType   dim  Type0 Type1 Type2 Type3 Type4 Type5")
    #     cc=0
    #     while cc<bs1:
    #         print(unpack(endian+'4s 8s i i i i i i i', buf[cc:cc+40]))
    #         cc+=40
    #     return(1)
    if block == 'IDTP':  # Particle type
        npf = open(filename, 'rb')
        if fmt != 0:
            bname, bsize = read_bhead(npf)
        bs1 = npf.read(4)  # size of header
        npart = np.zeros(6, dtype='int32')
        npart[:] = unpack(endian + 'i i i i i i', npf.read(4 * 6))
        npf.close()

        idtype = np.zeros(npart.sum(), dtype=np.int32)
        nn = 0
        for i, j in enumerate(npart):
            if j > 0:
                idtype[nn:nn + j] = i
                nn += j
        return(idtype)
    else:
        if fmt >= 0:
            npf = open(filename, 'rb')
            if fmt == 1:
                bname, bsize = read_bhead(npf)
            bs1 = npf.read(4)
            npart = np.zeros(6, dtype='int32')
            npart[:] = unpack(endian + 'i i i i i i', npf.read(4 * 6))

            if ptype == 0:
                pty = [0, npart[0]]
            elif ptype == 1:
                pty = [npart[0], npart[1]]
            elif ptype == 2:
                pty = [np.sum(npart[:2]), npart[2]]
            elif ptype == 3:
                pty = [np.sum(npart[:3]), npart[3]]
            elif ptype == 4:
                pty = [np.sum(npart[:4]), npart[4]]
            elif ptype == 5:
                pty = [np.sum(npart[:5]), npart[5]]
            else:
                pty = [0, np.sum(npart)]

            masstbl = np.zeros(6, dtype='float64')
            masstbl[:] = unpack(endian + 'd d d d d d', npf.read(8 * 6))
            if block == "MASS":
                idg0 = (npart > 0) & (masstbl <= 0)
                if (len(npart[idg0]) == 0) and (fullmass):  # No Mass block!
                    idg1 = (npart > 0) & (masstbl > 0)
                    if len(npart[idg1]) == 1:
                        return masstbl[idg1]
                    else:  # multi masstble
                        totmass = np.zeros(
                            np.sum(npart, dtype='int64'), dtype='float32')
                        countnm = 0
                        for i in np.arange(6):
                            if npart[i] > 0:
                                totmass[countnm:countnm +
                                        npart[i]] = masstbl[i]
                                countnm += npart[i]
                        return totmass
                if ptype is not None:
                    if masstbl[ptype] > 0:
                        return(masstbl[ptype])
            npf.close()

        npf = open(filename, 'rb')
        subdata = read_block(npf, block, endian, quiet,
                             longid, fmt, pty, rawdata)
        if subdata is not None:
            if block == "MASS":  # We fill the mass with the mass tbl value if needed
                npf.close()
                idg0 = (npart > 0) & (masstbl > 0)
                if (len(npart[idg0]) > 0) and (fullmass):
                    totmass = np.zeros(
                        np.sum(npart, dtype='int64'), dtype='float32')
                    bgc = 0
                    subc = 0
                    for k in np.arange(6):
                        if npart[k] > 0:
                            if(masstbl[k] > 0):
                                totmass[bgc:bgc + npart[k]
                                        ] = np.zeros(npart[k], dtype='float32') + masstbl[k]
                            else:
                                totmass[bgc:bgc + npart[k]
                                        ] = subdata[subc:subc + npart[k]]
                                subc += npart[k]
                            bgc += npart[k]
                    return totmass
                else:
                    if ptype is not None:
                        if (npart[ptype] == 0) or (masstbl[ptype] > 0):
                            print(
                                "This is can not be! npart[ptype] is ", npart[ptype], "masstbl[ptype] is ", masstbl[ptype])
                            print("I return 0")
                            return(0)
                        else:
                            startc = 0
                            endc = 0
                            for ii in range(ptype + 1):
                                if (npart[ii] > 0) and (masstbl[ii] <= 0):
                                    startc = endc
                                    endc += npart[ii]
                            return(subdata[startc:endc])
                    return subdata
            elif (block == "Zs  ") & (met == 'z'):
                if ptype == 0:
                    if masstbl[0] > 0:
                        mass = np.zeros(
                            npart[0], dtype=masstbl.dtype) + masstbl[0]
                    else:
                        mass = read_block(npf, "MASS", endian, 1, longid, fmt, [
                                          0, npart[0]], rawdata)
                    npf.close()
                    return np.sum(subdata[0:npart[0], 1:], axis=1) / (mass[0:npart[0]] - np.sum(subdata[0:npart[0], :], axis=1))
                elif ptype == 4:
                    im = read_block(npf, "iM  ", endian, 1,
                                    longid, fmt, pty, rawdata)
                    npf.close()
                    return np.sum(subdata[npart[0]:, 1:], axis=1) / (im - np.sum(subdata[npart[0]:, :], axis=1))
                else:
                    zs = np.zeros(npart[0] + npart[4], dtype=subdata.dtype)
                    if masstbl[0] > 0:
                        mass = np.zeros(
                            npart[0], dtype=masstbl.dtype) + masstbl[0]
                    else:
                        mass = read_block(npf, "MASS", endian, 1, longid, fmt, [
                                          0, npart[0]], rawdata)
                    zs[0:npart[0]] = np.sum(subdata[0:npart[0], 1:], axis=1) / (
                        mass[0:npart[0]] - np.sum(subdata[0:npart[0], :], axis=1))
                    mass = 0
                    im = read_block(npf, "iM  ", endian, 1,
                                    longid, fmt, pty, rawdata)
                    zs[npart[0]:] = np.sum(
                        subdata[npart[0]:, 1:], axis=1) / (im - np.sum(subdata[npart[0]:, :], axis=1))
                    im, subdata = 0, 0
                    npf.close()
                    return zs
            elif (block == "Z   ") and (ptype is not None):
                npf.close()
                if ptype == 0:
                    return subdata[:npart[0]]
                elif ptype == 4:
                    return subdata[-npart[4]:]
                else:
                    print("Don't support this particle type %i with Z" % ptype)
            else:
                npf.close()
                return subdata
        else:
            print("No such blocks!!! or Not add in this reading!!!", block)
            npf.close()
            return(0)

# Read Block


def read_block(npf, block, endian, quiet, longid, fmt, pty, rawdata):
    global nmets
    endf = fstat(npf.fileno()).st_size

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

            if (block == 'POS ') & (loopnum == 0):
                return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
            elif (block == 'VEL ') & (loopnum == 1):
                return read_bdata(npf, 3, np.dtype('float32'), endian, pty)
            elif (block == 'ID  ') & (loopnum == 2):
                if longid:
                    return read_bdata(npf, 1, np.dtype('uint64'), endian, pty)
                else:
                    return read_bdata(npf, 1, np.dtype('uint32'), endian, pty)
            # elif (block=='MASS') & (loopnum==3):
            #    return read_bdata(npf,1,np.dtype('float32'),endian)
            elif (block == 'U   ') & (loopnum == 3):
                return read_bdata(npf, 1, np.dtype('float32'), endian)
            elif loopnum > 4:
                return None
            loopnum += 1

        if not quiet:
            if fmt != 0:
                print(bname, bsize)
            else:
                print("Format 0, reading block ", block,  "skiping", bsize)

        #  For reading snapshot files###
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
        elif bname == block == 'Z   ':  # My own block, which saves the metallicity z
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'CLDX':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'TEMP':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'HSML':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        elif bname == block == 'PTYP':
            return read_bdata(npf, 1, np.dtype('int32'), endian)
        #  Internal energy###
        elif bname == block == 'U   ':
            return read_bdata(npf, 1, np.dtype('float32'), endian)
        #  INFO ####
        elif bname == block == 'INFO':  # This is print out not return in array
            bs1 = unpack(endian + 'i', npf.read(4))[0]
            buf = npf.read(bs1)
            print("Block   DataType   dim  Type0 Type1 Type2 Type3 Type4 Type5")
            cc = 0
            while cc < bs1:
                print(unpack(endian + '4s 8s i i i i i i i', buf[cc:cc + 40]))
                cc += 40
            return(1)

        #  For reading new subfind files###
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

        #  For reading my PIAO outputs ###
        elif (bname == block == 'GOFF') or (bname == block == 'GHED') or (bname == block == 'GSBL') \
                or (bname == block == 'GSBO') or (bname == block == 'SBLN') or (bname == block == 'SBOF'):
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
        arr = np.ndarray(shape=np.int32(
            bs1 / dt.itemsize), dtype=dt, buffer=buf)
    else:
        arr = np.ndarray(
            shape=(np.int32(bs1 / dt.itemsize / column), column), dtype=dt, buffer=buf)

    if endian == '=':
        return arr
    else:
        return arr.byteswap()
    # if unpack(endian+'i',npf.read(4))[0] != bs1:
    #    print("data part not consistent!!")
