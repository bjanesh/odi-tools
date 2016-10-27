import os
import sys
import odi_config as odi
from pyraf import iraf

def force_update_bpm(img, ota):
    image = odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])
    mask_name = odi.bppath+'reproj_mask_'+ota+'.'+str(img[16:])
    BPM = mask_name.replace('fits','pl')
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','BPM')
    iraf.imutil.hedit.setParam('value',BPM)
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')
    return

def make_bpms(img, ota):
    # for i in range(len(images)):
    #   for key in OTA_dictionary:
    # create string for mask fits name
    mask_name = odi.bppath+'mask_'+ota+'.'+str(img[16:])
    BPM = mask_name.replace('fits','pl')
    if not os.path.isfile(BPM):
        mask,gaps = odi.mask_ota(img,ota)
        hdu = odi.fits.PrimaryHDU(mask.astype(float))
        if not os.path.isfile(mask_name):
            hdu.writeto(mask_name,clobber=True)
        #if not os.path.isfile(mask_name.replace('fits','pl')):
            iraf.unlearn(iraf.imutil.imcopy)
            iraf.imutil.imcopy.setParam('input',mask_name)
            iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
            iraf.imutil.imcopy.setParam('verbose','no')
            iraf.imutil.imcopy(mode='h')
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',img+'['+ota+']')
    iraf.imutil.hedit.setParam('fields','BPM')
    iraf.imutil.hedit.setParam('value',BPM)
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')
    if os.path.isfile(mask_name):
        os.remove(mask_name)
    return

def get_dims(header_string):
    dims = header_string[header_string.index('['):
                         header_string.index(']')+1].split(',')
    dims1 = int(dims[0].strip('['))
    dims2 = int(dims[1].strip(']'))
    return dims1,dims2


def check_mask_dim(img,ota):
    mask_name = odi.bppath+'reproj_mask_'+ota+'.'+str(img[16:])
    reproj_img = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    bgsub_img = odi.bgsubpath+'bgsub_'+ota+'.'+str(img[16:])

    mask_header = iraf.imhead(mask_name,Stdout=1)[0]
    mask_dim1, mask_dim2 = get_dims(mask_header)

    reproj_header = iraf.imhead(reproj_img,Stdout=1)[0]
    reproj_dim1, reproj_dim2 = get_dims(reproj_header)

    bgsub_header = iraf.imhead(bgsub_img,Stdout=1)[0]
    bgsub_dim1, bgsub_dim2 = get_dims(bgsub_header)

    dim_status = ((mask_dim1 == reproj_dim1 == bgsub_dim1) &
                  (mask_dim2 == reproj_dim2 == bgsub_dim2))
    if dim_status == False:
        print 'mask and image dimensions to not match for: ', img, ota
        print 'mask:', mask_dim1,mask_dim2
        print 'reproj:', reproj_dim1,reproj_dim2
        print 'bgsub:', bgsub_dim1,bgsub_dim2
    else:
        print 'dimension test passed'
        print 'mask:', mask_dim1,mask_dim2
        print 'reproj:', reproj_dim1,reproj_dim2
        print 'bgsub:', bgsub_dim1,bgsub_dim2
    return dim_status

def main():
    pass

if __name__ == '__main__':
    main()
