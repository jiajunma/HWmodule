from fractions import Fraction as frac
from itertools import chain, combinations 
from bisect import bisect

def concat_strblocks(*args, sep=' '):                                      
    """                                                                    
    Combine the arguments into a single block of strings                   
    """                                                                    
    BLK = [ b.splitlines() for b in args]                                  
    v = max((0, *(len(b) for b in BLK)))                                   
    hBLK = [ max((0,*(len(r)  for r in b))) for b in BLK]                  
    res = []                                                               
    for i in range(v):                                                     
        l = [b[i]+' '*(hBLK[j]-len(b[i])) if i<len(b) else ' '*hBLK[j]     
                  for j,b in enumerate(BLK)]                               
        res.append(sep.join(l))                                            
    return '\n'.join(res)                                                  


'''
The RSK algorithm. 
'''
def RSK(p):
    '''Given a permutation p, spit out a pair of Young tableaux'''
    P = []; Q = []
    def insert(m, n=0):
        '''Insert m into P, then place n in Q at the same place'''
        for r in range(len(P)):
            if m >= P[r][-1]:
                P[r].append(m); Q[r].append(n)
                return
            c = bisect(P[r], m)
            P[r][c],m = m,P[r][c]
        P.append([m])
        Q.append([n])

    for i in range(len(p)):
        insert(p[i], i+1)
    return (P,Q)


def RSK2part(p):
    P,Q = RSK(p)
    return tuple(len(r) for r in P)


def str_part(part, s='*'):
    part = reg_part(part,reverse=True)
    return '\n'.join([s*r for r in part])

def reg_part(part, reverse=False):
    """
    Regularize the partition.
    It is a tuple of decresing sequence of positive integers
    """
    part = [x for x in part if x>0]
    part.sort(reverse=reverse)
    return tuple(part)

def part_transpose(part, reverse=False):            
    part=sorted([x for x in part if x>0])      
    if len(part) == 0:                         
        res = []                          
    else:   
        tpart = []                             
        for i in range(part[-1]):              
            ri  = len([x for x in part if x>i])
            tpart.append(ri)                   
        res = sorted(tpart, reverse=reverse)  
    return tuple(res)

def reg_W_repn(tau, reverse=False):                         
    """                                                    
    Regularize the W_n repn paramterized by bipartition    
    """                                                    
    ntau = (reg_part(tau[0], reverse=reverse),             
            reg_part(tau[1], reverse=reverse))             
    return ntau                                            

def symbol2repn(sym, reverse=False):    
    symL, symR = sym                                 
    tauL = tuple(lam-i for i, lam in enumerate(symL))
    tauR = tuple(lam-i for i, lam in enumerate(symR))
    return reg_W_repn((tauL, tauR), reverse=reverse) 

def repn2symbol(tau, rtype='C'):                                  
    tauL,tauR = reg_W_repn(tau, reverse=False)                    
    if rtype == 'C' or rtype == 'B':                              
        lL = max(len(tauL),len(tauR)+1)                           
        lR = lL-1                                                 
    elif rtype == 'D':                                            
        lR = lL = max(len(tauL),len(tauR))                        
    else:                                                         
        raise Exception('Wrong type',rtype)                       
    symL = tuple(i+lam for i, lam                                 
                 in enumerate(chain((0,)*(lL-len(tauL)),tauL)))   
    symR = tuple(i+lam for i, lam                                 
                 in enumerate(chain((0,)*(lR-len(tauR)),tauR)))   
    return (symL,symR)

def springer_part2symb(part, rtype='C'):
    part = sorted(part)
    if (rtype == 'C' and len(part)%2 == 1) or \
        (rtype == 'B' and len(part)%2 == 0) or \
        (rtype == 'D' and len(part)%2 == 1):
        part.insert(0, 0)
    pp = [lam+i for i, lam in enumerate(part)]
    pe, po = [],[]
    for lam in pp:
        if lam % 2 == 1:
            po.append(lam//2)
        else:
            pe.append(lam//2)
    return (po,pe)
    #tauL = tuple(xis - i for i, xis in enumerate(po))
    #tauR = tuple(eta - i for i, eta in enumerate(pe))
    #return (tauL,tauR)
    
    
'''
Not correct for D
'''
def cell_part2symb(part, rtype='C'):
    part = sorted(part)
    if  len(part)%2 == 0:
        part.insert(0, 0)
    pp = [lam+i for i, lam in enumerate(part)]
    pe, po = [],[]
    for lam in pp:
        if lam % 2 == 1:
            po.append(lam//2)
        else:
            pe.append(lam//2)
    if rtype == 'C':
        symb = (pe,po)
    else:
        symb = (pe, (0, *(mu+1 for mu in po)))
    return symb



def specialsymbol(sym):          
    return ssymbol2symbol(chain(sym[0],sym[1]))  

def ssymbol2symbol(ssym):                       
    """                                         
    From the set of element to special symbol   
    """                                         
    ssym = sorted(ssym)                         
    tauL,tauR = [],[]                           
    for i in range(len(ssym)):                  
        if i%2 ==0 :                            
            tauL.append(ssym[i])                
        else:                                   
            tauR.append(ssym[i])                
    return (tauL,tauR)                          
                                                

def springer_repn2part(tau, rtype = 'C'):                               
    tauL, tauR = tau                                                    
    if rtype == 'C':                                                    
        lL = max(len(tauL),len(tauR)+1)                                 
        lR = lL-1                                                       
        tauL = sorted([0]*(lL-len(tauL))+list(tauL))                    
        tauR = sorted([0]*(lR-len(tauR))+list(tauR))                    
        xis = [xi+2*i for i, xi in enumerate(tauL)]                     
        etas = [eta+2*i+1 for i, eta in enumerate(tauR)]                
        ssym = sorted(xis+etas)                                         
        ssymL, ssymR = ssymbol2symbol(ssym)                             
        """                                                             
        symbol = (xi_i + 2i; eta_i + 2i+1)                              
        """                                                             
        sxi = [lam-2*i for i, lam in enumerate(ssymL)]                  
        seta = [0] + [lam-2*i-1 for i,lam in enumerate(ssymR)]          
        """                                                             
        Compute the partition                                           
        """                                                             
        olams = [(lam+i)*2+1 for i, lam in enumerate(sxi)]              
        elams = [(lam+i)*2 for i, lam in enumerate(seta)]               
        part = [lam - i for i, lam in enumerate(sorted(olams+elams))]   
        return reg_part(part)                                           
    elif rtype == 'B' or rtype == 'D':                                  
        if rtype == 'B':                                                
            lL = max(len(tauL),len(tauR)+1)                             
            lR = lL-1                                                   
        else:                                                           
            lR = lL = max(len(tauL),len(tauR))                          
        tauL = sorted([0]*(lL-len(tauL))+list(tauL))                    
        tauR = sorted([0]*(lR-len(tauR))+list(tauR))                    
        xis = [xi+2*i for i, xi in enumerate(tauL)]                     
        etas = [eta+2*i for i, eta in enumerate(tauR)]                  
        ssym = sorted(xis+etas)                                         
        ssymL, ssymR = ssymbol2symbol(ssym)                             
        """                                                             
        symbol = (xi_i + 2i; eta_i + 2i+1)                              
        """                                                             
        sxi = [lam-2*i for i, lam in enumerate(ssymL)]                  
        seta = [lam-2*i for i,lam in enumerate(ssymR)]                  
        """                                                             
        Compute the partition                                           
        """                                                             
        olams = [(lam+i)*2+1 for i, lam in enumerate(sxi)]              
        elams = [(lam+i)*2 for i, lam in enumerate(seta)]               
        part = [lam - i for i, lam in enumerate(sorted(olams+elams))]   
        return reg_part(part)                                           
    
def HWM(HW):
    return (*HW, *(-wt for wt in HW[::-1]))

def MHW(HW):
    return (*(-wt for wt in HW[::-1]), *HW)
    

        
def infsumA(A):
    return sum( min(a,ap)  for i,a in enumerate(A) for ap in A[i+1:])

def infsumAB(A,B):
    return sum(min(a,b)  for a in A for b in B)

def msum(m):
    return sum( a*(a-1)//2 for a in range(m,0,-2))

def repn2fakedegree(tau,rtype):
    '''
    We follow Carter's book Section 11.4
    '''
    symU, symD = repn2symbol(tau,rtype)
    m = len(symD)
    res = 0
    if rtype in ('B','C'):
        res += 2*infsumA(symU)+2*infsumA(symD)
        res += sum(symD)
        res -= msum(2*m-1)
    else:
        res += 2*infsumA(symU)+2*infsumA(symD)
        res += min(sum(symU),sum(symD))
        res -= msum(2*m-2)
    return res

def repn2genericdegree(tau,rtype):
    '''
    We follow Carter's book Section 11.4
    '''
    symU, symD = repn2symbol(tau,rtype)
    m = len(symD)
    res = 0
    if rtype in ('B','C'):
        res += infsumA(symU)+infsumA(symD)
        res += infsumAB(symU,symD)
        res -= msum(2*m-1)
    else:
        res += infsumA(symU)+infsumA(symD)
        res += infsumAB(symU,symD)
        res -= msum(2*m-2)
    return res
    
def jindD2B(tau):
    '''
    j-induction from D_n to BC_n:
    compute the fakedegree of (tauL,tauR) and (tauR,tauL) of W_n repn. 
    Pick the one has minimal fakedegree. 
    '''
    tauL, tauR = tau
    gd1 = repn2fakedegree(tau,'C')
    gd2 = repn2fakedegree((tauR,tauL),'C')
    if gd1>=gd2:
        restau = (tauR,tauL)
    else:
        restau = tau
    assert(repn2fakedegree(restau,'C') == repn2fakedegree(tau,'D'))
    #print(f'j-ind D2B init bipart: {tau}')
    #print(f'result bipart: {restau}')   
    return restau

def jindA2B(part):
    '''
    j-induction from A_n to BC_n
    part ==> dual partition = (a_1,a_2, ..., a_k)
    '''
    tpart = part_transpose(part)
    ttauL, ttauR = tuple((a+1)//2 for a in tpart if a>0), tuple(a//2 for a in tpart if a>0)
    tauL, tauR = part_transpose(ttauL), part_transpose(ttauR)
    #print(f'A2B init part: \n {str_part(part)}')
    #strbpart = concat_strblocks(str_part(tauL),',',str_part(tauR))
    #print(f'result bipart: \n {strbpart}')
    return (tauL,tauR)

def jindBS2B(TAUS):
    tTAUSL = sorted([a for tauL,tauR in TAUS for a in part_transpose(tauL)])
    tTAUSR = sorted([b for tauL,tauR in TAUS for b in part_transpose(tauR)])
    tauL,tauR = part_transpose(tTAUSL), part_transpose(tTAUSR)
    
    #print(f'TAUS: {TAUS}')
    #print(f'res: ({tauL}, {tauR})')
    return (tauL,tauR)


def jindBC(subrepns):
    BS = []
    for tau, rtype in subrepns:
        if rtype == 'A':
            #print(jindA2B(tau))
            BS.append(jindA2B(tau))
        elif rtype == 'D':
            BS.append(jindD2B(tau))
        elif rtype in ('B','C'):
            BS.append(tau)
        else:
            raise ValueError(f"{tau},{rtype}")
    return jindBS2B(BS)


'''
Compute the orbit dimension following the foluma in CM. 
'''
def part2rowtuple(part):
    '''
    translate [row_1, row_2, ... row_n] to tuple (r_1, r_2, ... r_{row_1})
    where r_i is the number of rows of lenght i.
    '''
    spart = part_transpose(part,reverse=True)
    if not spart:
        return tuple()
    else:
        return (*(spart[i]-spart[i+1] for i in range(0,len(spart)-1)), spart[-1])

def codimO(part, rtype='A'):
    S = part_transpose(part,reverse=True)
    if not S:
        return 0
    res = sum(s*s for s in S)
    if rtype == 'A':
         res = res - 1
    else:
        RR = part2rowtuple(part)  
        r = sum(r for r in RR[::2])
        if rtype == 'C':
            res = (res + r)//2
        elif rtype in ('B','D'):
            res = (res - r)//2
        else:
            raise ValueError(f"{part}, {rtype}")
    return res

def dimO(part, rtype='A'):
    N = sum(part)
    cd = codimO(part,rtype)
    res = 0
    if rtype == 'A':
        res = N*N - cd
    elif rtype == 'C':
        res = (N*(N+1) // 2) -cd
    elif rtype in ('B','D'):
        res = (N*(N-1) // 2)- cd 
    else:
        raise ValueError
    return res


IWSYMB = {'C':'C', 'B':'C','D':'D'}
HIWSYMB = {'C':'D', 'B':'C', 'D':'D'}

def HW2AV(HW, rtype='C', report = 0):
    '''
    highest weight + rho ===> complex associated variety 
    '''
    assert(rtype in ('C','B','D'))
    
    HW = tuple(frac(wt) for wt in HW)
    # dict of integral weights
    IW = {frac('0'):[], frac('1/2'):[]}
    for a in HW:
        r = frac(a)%1
        nr = frac(-a)%1
        if r == frac('0') or r == frac('1/2'):
            IW[r].append(a)
        elif r in IW:
            IW[r][0].append(a)
        elif nr in IW:
            IW[nr][1].append(a)
        else:
            IW[r] = [[a],[]]

    subrepns = []
    #integer weight
    iw = IW.pop(frac('0'))
    extiw = (*iw,*(-wt for wt in iw[::-1]))
    if report>2:
        print(f'integer weight: {iw}')
    iwpart = RSK2part(extiw)
    #print(str_part(iwpart))
    iwsym = specialsymbol(cell_part2symb(iwpart,rtype = IWSYMB[rtype]))
    iwtau = symbol2repn(iwsym)
    subrepns.append((iwtau, IWSYMB[rtype]))
    
    #halfinteger weight
    hiw = IW.pop(frac('0.5'))
    if report>=2:
        print(f'half integer weight: {hiw}')        
    exthiw = (*hiw,*(-wt for wt in hiw[::-1]))
    hiwpart = RSK2part(exthiw)
    #print(str_part(iwpart))
    hiwsym = specialsymbol(cell_part2symb(hiwpart, rtype = HIWSYMB[rtype]))
    hiwtau = symbol2repn(hiwsym)
    subrepns.append((hiwtau, HIWSYMB[rtype]))
    
    # non integral or half-integral weight
    for r, WPM in IW.items():
        wt = (*WPM[0], *(-t for t in WPM[1][::-1]))
        atau = RSK2part(wt)
        if report>=2:
            print(f'weights (mod {r}): {wt}')
            #print(str_part(atau))
        subrepns.append((atau,'A'))
        
    # j-induction to the big group
    tau = jindBC(subrepns)
    if report >=2:
        print(f'sub-repn. : {subrepns}')
        print(f'tau: {tau}')
    orbit = springer_repn2part(tau,rtype=rtype)
    orbit = reg_part(orbit,reverse=True)
    if report>=1:
        print(f'Type {rtype}_{sum(orbit)//2}:\n AV_C(L({F2S(HW)})) \n  = {orbit} ')
        print(' AVC-dim=',dimO(orbit,rtype=rtype))
        print(str_part(orbit))
    return orbit
                
def S2F(s):
    '''
    convert string to tuple of fractions
    '''
    return tuple(frac(t) for t in s.split(','))

def F2S(f):
    '''
    convert tuple of fractions to string
    '''
    return ','.join(str(t) for t in f) 


