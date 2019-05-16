

def get_param(f):
    line = f.readline()
    i = line.index(" ")
    return line[:i]

def read_params():
    f = open('params.txt', 'r')
    
    s = get_param(f); tmax  = float(s)
    s = get_param(f); nt    = int(s)
    s = get_param(f); beta  = float(s)
    s = get_param(f); ntau  = int(s)
    s = get_param(f); norb  = int(s)
    s = get_param(f); order = int(s)
    
    f.close()

    return tmax, nt, beta, ntau, norb, order

if __name__=='__main__':
    read_params()
