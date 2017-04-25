import sys, warnings, urllib.request

def readFromUrl(Url=None):
    #This code opens the file from the given Url link using internet connection.
    
    if Url == '':
        warnings.warn('Empty string was given as the Url!')
        sys.exit(0) #Check if string is empty
    else:
        req = urllib.request.Request(Url)
        try:
            with urllib.request.urlopen(req) as response:
                html = response.read()
        except urllib.error.HTTPError as e:
            warnings.warn('\n{}\n{}'.format(e.code, e.read()))
            sys.exit(0)
    
    UTF = html.decode('utf-8')
    FILE = str(UTF)[:-1].split('\n')
    return(FILE)