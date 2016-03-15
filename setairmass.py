###########
# run setairmass

if not os.path.isfile('setairmass.done'):
    iraf.astutil.setairmass.setParam('images', "msc*fits")          # Input images
    iraf.astutil.setairmass.setParam('intype', "beginning")    # Input keyword time stamp
    iraf.astutil.setairmass.setParam('outtype', "effective")    # Output airmass time stamp\n
    iraf.astutil.setairmass.setParam('ra', "ra")           # Right acsension keyword (hours)
    iraf.astutil.setairmass.setParam('dec', "dec")          # Declination keyword (degrees)
    iraf.astutil.setairmass.setParam('equinox', "radeceq")        # Equinox keyword (years)
    iraf.astutil.setairmass.setParam('st', "st")           # Local siderial time keyword (hours)
    iraf.astutil.setairmass.setParam('ut', "time-obs")     # Universal time keyword (hours)
    iraf.astutil.setairmass.setParam('date', "date-obs")     # Observation date keyword
    iraf.astutil.setairmass.setParam('exposure', "exptime")      # Exposure time keyword (seconds)
    iraf.astutil.setairmass.setParam('airmass', "airmass")      # Airmass keyword (output)
    iraf.astutil.setairmass.setParam('utmiddle', "utmiddle")     # Mid-observation UT keyword (output)
    iraf.astutil.setairmass.setParam('scale', 750.)           # The atmospheric scale height\n
    iraf.astutil.setairmass.setParam('show', 'yes')            # Print the airmasses and mid-UT?
    iraf.astutil.setairmass.setParam('update', 'yes')            # Update the image header?
    iraf.astutil.setairmass.setParam('override', 'yes')            # Override previous assignments?
    iraf.astutil.setairmass()
    with open('setairmass.done','w+') as f1:
        print >> f1, True
else:
    print 'setairmass already done'
