import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
from scipy.io import mmread
import networkx.algorithms.isomorphism as iso
import math
import sys
from operator import itemgetter
import pandas as pd
from networkx.algorithms.approximation import clique
inFile = sys.argv[1]
k = sys.argv[2]
k=int(k)
Matrix = mmread(inFile)
B = Matrix.todense()
df = pd.DataFrame(B, range(1, B.shape[0] + 1), range(1, B.shape[1] + 1))
G = nx.from_pandas_adjacency(df)
print(G)
print(G.nodes())

H=G.copy()
D=[]
Sstar=[]
n=G.number_of_nodes()
if k<0 or k>n:
    print('Invalid value of k')

class SVD:
    def __init__(self,G,k,D):
        self.G=G
        self.k=k
        self.D=D

    def RR1(self):
        H1=G.copy()
        #cycle-4 subgraph to check isomorphism
        cycle4_graph = nx.Graph()
        cycle4_graph.add_nodes_from([7, 10])
        cycle4_graph.add_edges_from([(7,8),(8,9),(9,10),(10,7)])

        #cycle-5 subgraph to check isomorphism
        cycle5_graph = nx.Graph()
        cycle5_graph.add_nodes_from([17, 21])
        cycle5_graph.add_edges_from([(17,18),(18,19),(19,20),(20,21),(21,17)])

        #Two K2 subgraph to check isomorphism
        twok2_graph = nx.Graph()
        twok2_graph.add_nodes_from([17, 20])
        twok2_graph.add_edges_from([(17,18),(19,20)])


        D1=[]
        D2=[]
        D3=[]
        i=0                                                                                 
        j=0 
        l=0  
        #RR1 : to find forbidden structures                                                                 
                                                                                
                                                                                #Cycle-5 is found and delete the nodes from G

        GM2 = iso.GraphMatcher(H1, cycle5_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM2.subgraph_is_isomorphic():
            print('Cycle with 5 vertices found')
            dic2=GM2.mapping
            #print(dic1)
            D2.append(list(dic2.keys()))
            #print(D2)
            H1.remove_nodes_from(D2[j])
            j=j+1
            #print(G.nodes())


                                                                                #Cycle-4 is found and delete the nodes from G
        GM1 = iso.GraphMatcher(H1, cycle4_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM1.subgraph_is_isomorphic():
            print('Cycle with 4 vertices found')
            dic1=GM1.mapping
            #print(dic1)
            D1.append(list(dic1.keys()))
            #print(D1)
            H1.remove_nodes_from(D1[i])
            i=i+1
            #print(G.nodes())


                                                                                #2K2 is found and delete the nodes from G    
        GM3 = iso.GraphMatcher(H1, twok2_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM3.subgraph_is_isomorphic():
            print('2K2 is found')
            dic3=GM3.mapping
            #print(dic3)
            D3.append(list(dic3.keys()))
            #print(D3)
            H1.remove_nodes_from(D3[l])
            l=l+1

        D4=D1+D2+D3
        self.D = [item for sublist in D4 for item in sublist]
        self.D.sort()
        print('D contain nodes:',self.D)     
        print('G-D contain nodes:',H1.nodes())
        print('Reduction Rule 1 completed!')
        print('\n')                                                                    #if |D|>5k it is a NO instance
        if (len(D)>(5*self.k)):
            print('NO Instance because |D| > 5k')
            sys.exit()
        else:
            return (self.G,self.k,self.D)

    
    def RR2(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        Cstar=[]
        Cstar.append(V[0])                                                            # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                            #Istar is I* ===> C* \disjoint_union\ I* is a split graph

        #RR2 begins
        
        
        Hi=[]
        Hc=[]
        for v in V:
            if (len([n for n in Cstar if n not in self.G.neighbors(v)]) >= (self.k+2)):              # Hi={x \in V: |C*\N(v)|>=k+2}
                Hi.append(v) 
            if (len([n for n in Istar if n in self.G.neighbors(v)]) >= (self.k+2)):
                Hc.append(v)
        print('Hi=',Hi)
        print('Hc=',Hc)

        # RR2- if there is a vertex in Hi \intersection\ Hc, delete v and reduce k by 1

        vinter=[n for n in Hi if n in Hc]
        if len(vinter)>0:
            Sstar.extend(vinter)
            self.k=self.k-len(vinter)
            self.G.remove_nodes_from(vinter)
            
        print(self.G.nodes)    

        print("Reduction Rule 2 completed! ")
        print('\n')
        return (self.G,self.k,self.D)
        
    
    def RR3(self):
                                                              #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=list(self.G.nodes())
        E=self.G.edges()
        Cstar=[]
        Cstar.append(V[0])
        #print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        #print('I*=',Istar)  
        Hi=[]
        Hc=[]
        for v in V:
            if (len([n for n in Cstar if n not in self.G.neighbors(v)]) >= (self.k+2)):              # Hi={x \in V: |C*\N(v)|>=k+2}
                Hi.append(v) 
            if (len([n for n in Istar if n in self.G.neighbors(v)]) >= (self.k+2)):
                Hc.append(v)
        print('Hi=',Hi)
        print('Hc=',Hc)
        C1=[n for n in Hc if n in Cstar]                    #C1 = Hc ∩ C*
        I1=[n for n in Hi if n in Istar]                    #I1 = Hi ∩ I*
        C0=[n for n in Hc if n in self.D]                        #C0 = Hc ∩ D
        I0=[n for n in Hi if n in self.D]                        #I1 = Hi ∩ D
        C1star=[n for n in Cstar if n not in C1]            #C1star = Cstar \ C1 
        I1star=[n for n in Istar if n not in I1]            #I1star = Istar \ I1
        Y1=[]
        Y2=[]
        DminusI0=[n for n in self.D if n not in I0]
        DminusI0UI1star=DminusI0 + I1star
        for v in C1star:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y1.append(v)
        for v in C1:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y2.append(v)
        C1r=[n for n in C1star if n not in Y1]
        C2r=[n for n in C1 if n not in Y2]
        X1=[]
        X2=[]
        DminusC0=[n for n in self.D if n not in C0]
        DminusC0UC1star=DminusC0 + C1star
        for v in I1star:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X1.append(v)
        for v in I1:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X2.append(v)
        I1r=[n for n in I1star if n not in X1]
        I2r=[n for n in I1 if n not in X2]

        #RR3 - If |C1r|>k+2, delete all edges between C1r and I1 U I0 and delete all but k+2 vertices of C1r

        I1UI0=I1+I0

        if len(C1r) > (self.k)+2:
            for v in C1r:
                vnbr=[n for n in self.G.neighbors(v) if n in I1UI0]
                if len(vnbr)>0:
                    for u in vnbr:
                        self.G.remove_edge(v,u)
            for i in C1r:
                if len(C1r)==(self.k)+2:
                    break
                C1r.remove(i)
        print('Reduction Rule 3 completed!')
        print('\n')
        return (self.G,self.k,self.D)


    def RR4(self):
                                                              #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=list(self.G.nodes())
        E=self.G.edges()
        Cstar=[]
        Cstar.append(V[0])
        #print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        #print('I*=',Istar)  
        Hi=[]
        Hc=[]
        for v in V:
            if (len([n for n in Cstar if n not in self.G.neighbors(v)]) >= (self.k+2)):              # Hi={x \in V: |C*\N(v)|>=k+2}
                Hi.append(v) 
            if (len([n for n in Istar if n in self.G.neighbors(v)]) >= (self.k+2)):
                Hc.append(v)
        print('Hi=',Hi)
        print('Hc=',Hc)
        C1=[n for n in Hc if n in Cstar]                    #C1 = Hc ∩ C*
        I1=[n for n in Hi if n in Istar]                    #I1 = Hi ∩ I*

        C0=[n for n in Hc if n in self.D]                        #C0 = Hc ∩ D
        I0=[n for n in Hi if n in self.D]                        #I1 = Hi ∩ D

        C1star=[n for n in Cstar if n not in C1]            #C1star = Cstar \ C1 
        I1star=[n for n in Istar if n not in I1]            #I1star = Istar \ I1

        Y1=[]
        Y2=[]

        DminusI0=[n for n in self.D if n not in I0]
        DminusI0UI1star=DminusI0 + I1star

        for v in C1star:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y1.append(v)

        for v in C1:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y2.append(v)

        C1r=[n for n in C1star if n not in Y1]
        C2r=[n for n in C1 if n not in Y2]


        X1=[]
        X2=[]

        DminusC0=[n for n in self.D if n not in C0]
        DminusC0UC1star=DminusC0 + C1star

        for v in I1star:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X1.append(v)

        for v in I1:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X2.append(v)

        I1r=[n for n in I1star if n not in X1]
        I2r=[n for n in I1 if n not in X2]
        I1UI0=I1+I0
        #RR4 - If |C2r| > k + 2, then delete all edges between C2r and I1 ∪ I0 and delete all but k + 2 vertices of C2r

        if len(C2r) > k+2:
            for v in C2r:
                vnbr=[n for n in G.neighbors(v) if n in I1UI0]
                if len(vnbr)>0:
                    for u in vnbr:
                        G.remove_edge(v,u)
            for i in C2r:
                if len(C2r)==k+2:
                    break
                C2r.remove(i)
        print('Reduction Rule 4 completed!')
        print('\n')
        return (self.G,self.k,self.D)


    def RR5(self):
                                                              #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=list(self.G.nodes())
        E=self.G.edges()
        Cstar=[]
        Cstar.append(V[0])
        #print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        #print('I*=',Istar)  
        Hi=[]
        Hc=[]
        for v in V:
            if (len([n for n in Cstar if n not in self.G.neighbors(v)]) >= (self.k+2)):              # Hi={x \in V: |C*\N(v)|>=k+2}
                Hi.append(v) 
            if (len([n for n in Istar if n in self.G.neighbors(v)]) >= (self.k+2)):
                Hc.append(v)
        print('Hi=',Hi)
        print('Hc=',Hc)
        C1=[n for n in Hc if n in Cstar]                    #C1 = Hc ∩ C*
        I1=[n for n in Hi if n in Istar]                    #I1 = Hi ∩ I*

        C0=[n for n in Hc if n in self.D]                        #C0 = Hc ∩ D
        I0=[n for n in Hi if n in self.D]                        #I1 = Hi ∩ D

        C1star=[n for n in Cstar if n not in C1]            #C1star = Cstar \ C1 
        I1star=[n for n in Istar if n not in I1]            #I1star = Istar \ I1

        Y1=[]
        Y2=[]

        DminusI0=[n for n in self.D if n not in I0]
        DminusI0UI1star=DminusI0 + I1star

        for v in C1star:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y1.append(v)

        for v in C1:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y2.append(v)

        C1r=[n for n in C1star if n not in Y1]
        C2r=[n for n in C1 if n not in Y2]


        X1=[]
        X2=[]

        DminusC0=[n for n in self.D if n not in C0]
        DminusC0UC1star=DminusC0 + C1star

        for v in I1star:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X1.append(v)

        for v in I1:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X2.append(v)

        I1r=[n for n in I1star if n not in X1]
        I2r=[n for n in I1 if n not in X2]

        C1UC0=C1+C0
        # #RR5 -  If |I1r| > k + 2, then add all edges between I1r and C1 U C0 and delete all but k+2 vertices of I1r
        if len(I1r) > k+2:
            for v in I1r:
                for u in C1UC0:
                    G.add_edge(v,u)
            for i in I1r:
                if len(I1r)==k+2:
                    break
                I1r.remove(i)
        print('Reduction Rule 5 Completed!')
        print('\n')
        return (self.G,self.k,self.D)



    def RR6(self):
                                                              #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=list(self.G.nodes())
        E=self.G.edges()
        Cstar=[]
        Cstar.append(V[0])                                                                # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar] 
        Hi=[]
        Hc=[]
        for v in V:
            if (len([n for n in Cstar if n not in self.G.neighbors(v)]) >= (self.k+2)):              # Hi={x \in V: |C*\N(v)|>=k+2}
                Hi.append(v) 
            if (len([n for n in Istar if n in self.G.neighbors(v)]) >= (self.k+2)):
                Hc.append(v)
        print('Hi=',Hi)
        print('Hc=',Hc)
        C1=[n for n in Hc if n in Cstar]                    #C1 = Hc ∩ C*
        I1=[n for n in Hi if n in Istar]                    #I1 = Hi ∩ I*

        C0=[n for n in Hc if n in self.D]                        #C0 = Hc ∩ D
        I0=[n for n in Hi if n in self.D]                        #I1 = Hi ∩ D

        C1star=[n for n in Cstar if n not in C1]            #C1star = Cstar \ C1 
        I1star=[n for n in Istar if n not in I1]            #I1star = Istar \ I1

        Y1=[]
        Y2=[]

        DminusI0=[n for n in self.D if n not in I0]
        DminusI0UI1star=DminusI0 + I1star
        for v in C1star:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y1.append(v)

        for v in C1:
            if len([n for n in DminusI0UI1star if n not in self.G.neighbors(v)]) > 0:
                Y2.append(v)

        C1r=[n for n in C1star if n not in Y1]
        C2r=[n for n in C1 if n not in Y2]


        X1=[]
        X2=[]

        DminusC0=[n for n in self.D if n not in C0]
        DminusC0UC1star=DminusC0 + C1star

        for v in I1star:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X1.append(v)

        for v in I1:
            if len([n for n in DminusC0UC1star if n in self.G.neighbors(v)]) > 0:
                X2.append(v)

        I1r=[n for n in I1star if n not in X1]
        I2r=[n for n in I1 if n not in X2]

        C1UC0=C1+C0
        #RR6 -  If |I2r| > k + 2, then add all edges between I2r and C1 U C0 and delete all but k+2 vertices of I2r
        if len(I2r) > k+2:
            for v in I2r:
                for u in C1UC0:
                    G.add_edge(v,u)
            for i in I2r:
                if len(I2r)==k+2:
                    break
                I2r.remove(i)
        
        print('Reduction Rule 6 Completed!')
        print('\n')
        return (self.G,self.k,self.D)


# #When none of the reduction rules apply, |C1*| ≤ 2k + 2 or |I1*| ≤ 2k + 2.

# #When none of the reduction rules apply, the number of vertices in C1* ∪ I1* is O(k^2).

# #When none of the reduction rules apply, the sets C* and I* contain O(k^3) vertices.



    def RR7(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        C=[]
        C.append(V[0])                                                                # Cstar is C*

        I=[n for n in self.G.nodes() if n not in C]
        S=self.D   
        
        #RR7 begins::If there is v ∈ SY , then delete v and decrease k by 1.
        SY=[]
        
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(nbrI) >= (self.k)+2 and len(non_nbrC) >=(self.k)+2 :
                SY.append(v) 
        if len(SY)>0:
            Sstar.extend(SY)
            self.k=self.k-len(SY)
            self.G.remove_nodes_from(SY)
            
        print('Reduction Rule 7 completed!')
        print('\n')
        return (self.G,self.k,self.D)



    def RR8(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        C=[]
        C.append(V[0])                                                                 # Cstar is C*

        I=[n for n in self.G.nodes() if n not in C] 
        SY=[]
        
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(nbrI) >= (self.k)+2 and len(non_nbrC) >=(self.k)+2 :
                SY.append(v)
        SC=[]
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            
            if len(nbrI) >= ( self.k) +2 :
                SC.append(v) 
        SI=[]
        for v in self.D:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= (self.k) +2 :
                SI.append(v) 
        SCUSI=SC+SI
        SZ=[n for n in self.D if n not in SCUSI]
        Cfix=[]
        for v in C:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= k+2 :
                Cfix.append(v) 
        CZ=[n for n in C if n not in Cfix]
        Ifix=[]
        for v in I:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= k+2 :
                Ifix.append(v) 
        IZ=[n for n in I if n not in Ifix]
        SCUSZ=SC + SZ 
        SCUSZUIZ=SCUSZ+IZ
        CY=[]
        for v in C:
            no_non_nbrs=[n for n in SCUSZUIZ if n not in self.G.neighbors(v)]
            if len(no_non_nbrs)==0:
                CY.append(v)
        CYCfix=[n for n in CY if n in Cfix]
        IfixUSI= Ifix + SI

        #RR8 begins:: If |CY ∩ Cfix|>k+2, then delete all edges between CY ∩ Cfix and (Ifix ∪ SI) and delete all but k+2 vertices of CY ∩ Cfix.
        if len(CYCfix) > k+2:
            for v in CYCfix:
                vnbr=[n for n in self.G.neighbors(v) if n in IfixUSI]
                if len(vnbr)>0:
                    for u in vnbr:
                        self.G.remove_edge(v,u)
            for i in CYCfix:
                if len(CYCfix)==(self.k)+2:
                    break
                CYCfix.remove(i)  
        print('Reduction Rule 8 completed!')
        print('\n')
        return (self.G,self.k,self.D)




    def RR9(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        C=[]
        C.append(V[0])                                                                # Cstar is C*

        I=[n for n in self.G.nodes() if n not in C] 
        SY=[]
        S=self.D
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(nbrI) >= (self.k)+2 and len(non_nbrC) >=(self.k)+2 :
                SY.append(v)
        SC=[]
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= ( self.k) +2 :
                SC.append(v) 
        SI=[]
        for v in self.D:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= (self.k) +2 :
                SI.append(v) 
        SCUSI=SC+SI
        SZ=[n for n in S if n not in SCUSI]
        Cfix=[]
        for v in C:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= k+2 :
                Cfix.append(v) 
        CZ=[n for n in C if n not in Cfix]
        Ifix=[]
        for v in I:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= k+2 :
                Ifix.append(v) 
        IZ=[n for n in I if n not in Ifix]
        SCUSZ=SC + SZ 
        SCUSZUIZ=SCUSZ+IZ
        CY=[]
        for v in C:
            no_non_nbrs=[n for n in SCUSZUIZ if n not in self.G.neighbors(v)]
            if len(no_non_nbrs)==0:
                CY.append(v)
        SIUSZ=SI + SZ 
        SIUSZUCZ=SIUSZ+CZ
        IY=[]
        for v in I:
            nbrs=[n for n in SIUSZUCZ if n in self.G.neighbors(v)]
            if len(nbrs)==0:
                IY.append(v)
        CYCfix=[n for n in CY if n in Cfix]
        IfixUSI= Ifix + SI

        #RR9 begins:: If |CY ∩ CZ |>k+2, then delete all edges between CY ∩ CZ and (Ifix ∪ SI) and delete all but k + 2 vertices of CY ∩ CZ .
        CYCZ=[n for n in CY if n in CZ]
        if len(CYCZ) > (self.k)+2:
            for v in CYCZ:
                vnbr=[n for n in self.G.neighbors(v) if n in IfixUSI]
                if len(vnbr)>0:
                    for u in vnbr:
                        self.G.remove_edge(v,u)
            for i in CYCZ:
                if len(CYCZ)==(self.k)+2:
                    break
                CYCZ.remove(i)
        print('Reduction Rule 9 completed!')
        print('\n')
        return (self.G,self.k,self.D)



    def RR10(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        C=[]
        C.append(V[0])                                                               # Cstar is C*

        I=[n for n in self.G.nodes() if n not in C]
        SY=[]
        S=self.D
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(nbrI) >= (self.k)+2 and len(non_nbrC) >=(self.k)+2 :
                SY.append(v)
        SC=[]
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= ( self.k) +2 :
                SC.append(v) 
        SI=[]
        for v in self.D:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= (self.k) +2 :
                SI.append(v) 
        SCUSI=SC+SI
        SZ=[n for n in S if n not in SCUSI]
        Cfix=[]
        for v in C:
            nbrI = [n for n in self.G.neighbors(v) if n in I]
            if len(nbrI) >= k+2 :
                Cfix.append(v) 
        CZ=[n for n in C if n not in Cfix]
        Ifix=[]
        for v in I:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= k+2 :
                Ifix.append(v) 
        IZ=[n for n in I if n not in Ifix]
        SCUSZ=SC + SZ 
        SCUSZUIZ=SCUSZ+IZ
        CY=[]
        for v in C:
            no_non_nbrs=[n for n in SCUSZUIZ if n not in self.G.neighbors(v)]
            if len(no_non_nbrs)==0:
                CY.append(v)
        SIUSZ=SI + SZ 
        SIUSZUCZ=SIUSZ+CZ
        IY=[]
        for v in I:
            nbrs=[n for n in SIUSZUCZ if n in self.G.neighbors(v)]
            if len(nbrs)==0:
                IY.append(v)
        CYCfix=[n for n in CY if n in Cfix]
        IfixUSI= Ifix + SI

        #RR10 :: If |IY ∩ Ifix| > k + 2, then add all edges between IY ∩ Ifix and (Cfix ∪ SC ) and delete all but k + 2 vertices of IY ∩ Ifix.

        IYIfix=[n for n in IY if n in Ifix]
        CfixUSC= Cfix + SC
        if len(IYIfix) > (self.k)+2:
            for v in IYIfix:
                for u in CfixUSC:
                    self.G.add_edge(v,u)
            for i in IYIfix:
                if len(IYIfix)==(self.k)+2:
                    break
                IYIfix.remove(i)
        print('Reduction Rule 10 completed!')
        print('\n')
        return (self.G,self.k,self.D)



    def RR11(self):
        V=list(self.G.nodes())
        E=self.G.edges()
        C=[]
        C.append(V[0])                                                                 # Cstar is C*

        I=[n for n in self.G.nodes() if n not in C]
        SY=[]
        S=self.D
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(nbrI) >= (self.k)+2 and len(non_nbrC) >=(self.k)+2 :
                SY.append(v)
        SC=[]
        for v in self.D:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= ( self.k) +2 :
                SC.append(v) 
        SI=[]
        for v in self.D:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= (self.k) +2 :
                SI.append(v) 
        SCUSI=SC+SI
        SZ=[n for n in S if n not in SCUSI]
        Cfix=[]
        for v in C:
            try:
                nbrI = [n for n in I if n in self.G.neighbors(v)]
            except:
                nbrI=[]
            if len(nbrI) >= k+2 :
                Cfix.append(v) 
        CZ=[n for n in C if n not in Cfix]
        Ifix=[]
        for v in I:
            try:
                non_nbr_v = [n for n in self.G.nodes() if n not in self.G.neighbors(v)]
            except:
                non_nbr_v=self.G.nodes()
            non_nbrC = [n for n in non_nbr_v if n in C]
            if len(non_nbrC) >= k+2 :
                Ifix.append(v) 
        IZ=[n for n in I if n not in Ifix]
        SCUSZ=SC + SZ 
        SCUSZUIZ=SCUSZ+IZ
        CY=[]
        for v in C:
            no_non_nbrs=[n for n in SCUSZUIZ if n not in self.G.neighbors(v)]
            if len(no_non_nbrs)==0:
                CY.append(v)
        SIUSZ=SI + SZ 
        SIUSZUCZ=SIUSZ+CZ
        IY=[]
        for v in I:
            nbrs=[n for n in SIUSZUCZ if n in self.G.neighbors(v)]
            if len(nbrs)==0:
                IY.append(v)
        CYCfix=[n for n in CY if n in Cfix]
        IfixUSI= Ifix + SI

        #RR11 :: If |IY ∩ I Z | > k + 2, then add all edges between IY ∩ I Z and (Cfix ∪ SC ) and delete all but k + 2 vertices of IY ∩ I Z .

        CfixUSC=Cfix +SC
        IYIZ=[n for n in IY if n in IZ]
        if len(IYIZ) > k+2:
            for v in IYIZ:
                for u in CfixUSC:
                    G.add_edge(v,u)
            for i in IYIZ:
                if len(IYIZ)==k+2:
                    break
                IYIZ.remove(i)
        print('Reduction Rule 11 completed!')
        print('\n')
        return (self.G,self.k,self.D)







#stage 1
obj=SVD(G,k,D)
G1,k,D=obj.RR1()

#stage 2
obj=SVD(G1,k,D)
G2,k,D=obj.RR2()

obj=SVD(G2,k,D)
G3,k,D=obj.RR1()

#stage 3
obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR2()

obj=SVD(G5,k,D)
G6,k,D=obj.RR1()

#stage 4
obj=SVD(G6,k,D)
G7,k,D=obj.RR4()

obj=SVD(G7,k,D)
G8,k,D=obj.RR1()

obj=SVD(G8,k,D)
G9,k,D=obj.RR2()

obj=SVD(G9,k,D)
G10,k,D=obj.RR3()

#stage 5
obj=SVD(G10,k,D)
G1,k,D=obj.RR5()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()


#stage 6
obj=SVD(G5,k,D)
G1,k,D=obj.RR6()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G4,k,D)
G6,k,D=obj.RR3()

obj=SVD(G6,k,D)
G7,k,D=obj.RR4()

obj=SVD(G7,k,D)
G8,k,D=obj.RR5()


#stage 7
obj=SVD(G8,k,D)
G1,k,D=obj.RR7()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()

obj=SVD(G5,k,D)
G6,k,D=obj.RR5()

obj=SVD(G6,k,D)
G7,k,D=obj.RR6()


#stage 8
obj=SVD(G7,k,D)
G1,k,D=obj.RR8()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()

obj=SVD(G5,k,D)
G6,k,D=obj.RR5()

obj=SVD(G6,k,D)
G7,k,D=obj.RR6()

obj=SVD(G7,k,D)
G8,k,D=obj.RR7()

#stage 9
obj=SVD(G8,k,D)
G1,k,D=obj.RR9()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()

obj=SVD(G5,k,D)
G6,k,D=obj.RR5()

obj=SVD(G6,k,D)
G7,k,D=obj.RR6()

obj=SVD(G7,k,D)
G8,k,D=obj.RR7()

obj=SVD(G8,k,D)
G9,k,D=obj.RR8()


#stage 10
obj=SVD(G9,k,D)
G1,k,D=obj.RR10()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()

obj=SVD(G5,k,D)
G6,k,D=obj.RR5()

obj=SVD(G6,k,D)
G7,k,D=obj.RR6()

obj=SVD(G7,k,D)
G8,k,D=obj.RR7()

obj=SVD(G8,k,D)
G9,k,D=obj.RR8()

obj=SVD(G9,k,D)
G10,k,D=obj.RR9()


#stage 11
obj=SVD(G10,k,D)
G1,k,D=obj.RR11()

obj=SVD(G1,k,D)
G2,k,D=obj.RR1()

obj=SVD(G2,k,D)
G3,k,D=obj.RR2()

obj=SVD(G3,k,D)
G4,k,D=obj.RR3()

obj=SVD(G4,k,D)
G5,k,D=obj.RR4()

obj=SVD(G5,k,D)
G6,k,D=obj.RR5()

obj=SVD(G6,k,D)
G7,k,D=obj.RR6()

obj=SVD(G7,k,D)
G8,k,D=obj.RR7()

obj=SVD(G8,k,D)
G9,k,D=obj.RR8()

obj=SVD(G9,k,D)
G10,k,D=obj.RR9()

obj=SVD(G10,k,D)
G11,k,D=obj.RR10()


def rho(Gn,v):
    nonnbrv=[n for n in Gn.nodes() if n not in Gn.neighbors(v)]
    nonnbrvCU=[n for n in nonnbrv if n in CU]
    return len(nonnbrvCU)

print('Graph after applying reduction rules:',G11.nodes())
print('\n')

Capp=nx.Graph()
maxi=max(H.nodes())
Capp.add_nodes_from(range(maxi+1,maxi+1+k+2))
fromnodes=range(maxi+1,maxi+1+k+2)
tonodes=range(maxi+1,maxi+1+k+2)
for x in fromnodes:
    for y in tonodes:
        if x!=y:
            Capp.add_edge(x,y)
Iapp=nx.Graph()
Iapp.add_nodes_from(range(maxi+k+3,maxi+k+3+k+2))


V=list(G11.nodes())
C=[]
C.append(V[0])
print('C=',C)                                                                 # Cstar is C*

I=[n for n in G11.nodes() if n not in C]                                     
print('I=',I)  
S=D
SC=[]
for v in D:
    try:
        nbrI = [n for n in I if n in G11.neighbors(v)]
    except:
        nbrI=[]
    if len(nbrI) >= ( k) +2 :
        SC.append(v) 
SI=[]
for v in D:
    try:
        non_nbr_v = [n for n in G11.nodes() if n not in G11.neighbors(v)]
    except:
        non_nbr_v=G11.nodes()
    non_nbrC = [n for n in non_nbr_v if n in C]
    if len(non_nbrC) >= (k) +2 :
        SI.append(v) 
SCUSI=SC+SI
SZ=[n for n in S if n not in SCUSI]
Cfix=[]
for v in C:
    try:
        nbrI = [n for n in I if n in G11.neighbors(v)]
    except:
        nbrI=[]
    if len(nbrI) >= k+2 :
        Cfix.append(v) 
CZ=[n for n in C if n not in Cfix]
Ifix=[]
for v in I:
    try:
        non_nbr_v = [n for n in G11.nodes() if n not in G11.neighbors(v)]
    except:
        non_nbr_v=G11.nodes()
    non_nbrC = [n for n in non_nbr_v if n in C]
    if len(non_nbrC) >= k+2 :
        Ifix.append(v) 
IZ=[n for n in I if n not in Ifix]

CUSC=C+SC
CUSCUSZ=CUSC+SZ
CUSCUSZUIZ=CUSCUSZ+IZ
Gapp=nx.Graph()

G1=nx.compose(G11,Capp)

Gapp=nx.compose(G1,Iapp)

for u in Capp.nodes():
    for v in CUSCUSZUIZ:
        if u in Gapp.nodes() and v in Gapp.nodes():
            Gapp.add_edge(u,v)

CfixUSC=Cfix+SC
CfixUSCUCapp=CfixUSC+list(Capp.nodes())

for u in Iapp.nodes():
    for v in CfixUSCUCapp:
        if u in Gapp.nodes() and v in Gapp.nodes():
            Gapp.add_edge(u,v)



print('Graph after appending a clique of k+2 vertices and an independent set of k+2 vertices:',Gapp.nodes(),len(Gapp))
print('\n')               #Now size of Gapp will be |V(G)|+2(k+2)

#M = Capp ∪ Iapp ∪ S ∪ I ∪ CZ

CappUIapp=list(Capp.nodes())+list(Iapp.nodes())
CappUIappUS=CappUIapp +S
CappUIappUSUI=CappUIappUS +I
M=CappUIappUSUI+CZ
SCUSZ=SC+SZ
Grho=Gapp.copy()
for s in SCUSZ:
    nonnbrs=[n for n in G11.nodes() if n not in G11.neighbors(s)]
    nonnbrsinC=[n for n in nonnbrs if n in C]
    M.append(nonnbrsinC)
CU=[n for n in Cfix if n not in M]

for v in CU:
    Grho.remove_node(v)


print('Graph after removing CU:',Gapp.nodes,len(Gapp))                 
print('\n')

Q=nx.Graph()
maxi=max(Gapp.nodes())
Q.add_nodes_from(range(maxi+1,maxi+1+k+1))
fromnodes=range(maxi+1,maxi+1+k+1)
tonodes=range(maxi+1,maxi+1+k+1)
for x in fromnodes:
    for y in tonodes:
        if x!=y:
            Q.add_edge(x,y)

Gf=nx.compose(Grho,Q)
#Capp ∪ Iapp ∪ (C \ CU ) ∪ SC ∪ S Z
CappUIapp=list(Capp.nodes())+list(Iapp.nodes())
CminusCU=[n for n in C if n not in CU]
CappUIappUCminusCU=CappUIapp + CminusCU
CappUIappUCminusCUUSCUSZ=CappUIappUCminusCU+SCUSZ
if len(CappUIappUCminusCUUSCUSZ)!=0:
    for q in Q:
        for v in CappUIappUCminusCUUSCUSZ:
            if q in Gf and v in Gf:
                Gf.add_edge(q,v)
Qlist=list(Q.nodes)
for v in IZ:
    ilist=[]
    for n in range(k+1):
        ilist.append(n)
    rhov=rho(G11,v)
    ilist.remove(rhov)
    for i in ilist:
        Gf.add_edge(v,Qlist[i])
print('Final graph after appending Q: ',Gf.nodes())   

#Algorithm that generate family P
nodes=list(Gf.nodes())


#print("li=",li)

def Generator(Gf,d,VC0 , VI0 , A,ans):
    n=Gf.number_of_nodes()
   
    VC0UA=VC0+A
    res=[]
    for i in VC0UA:
        if i not in res:
            res.append(i)
    VC0UA=res
    res=[]
    VI0UA=VI0+A
    for i in VI0UA:
        if i not in res:
            res.append(i)
    VI0UA=res
    ans.append([VC0UA , VI0 ])
    ans.append([VC0 , VI0UA])
    if (d < (2*(math.floor(math.log(n)+ 1)))):
        for v in A:
            Aminusngv=[i for i in A if i not in (Gf.neighbors(v))]
            if v in Aminusngv:
                Aminusngv.remove(v)
            Ainngv=[j for j in A if j in Gf.neighbors(v)]
            if v not in VC0:
                VC0.append(v)
            VI0UAminusngv=list(set(VI0) | set(Aminusngv))
            VC0UAinngv=list(set(VC0) | set(Ainngv))
            return Generator(Gf,d + 1,VC0, VI0UAminusngv ,Ainngv,ans)
            VC0UAinngv=list(set(VC0) | set(Ainngv))
            if v not in VI0:
                VI0.append(v)
            return Generator(Gf,d + 1,VC0UAinngv, VI0 ,Aminusngv,ans)
    else:
        return ans  

"""

def VertexCover(G,k,S):
    if G.edges()==[]:
        return True,S
    if k <=0:
        return False,S
    dlist=[]
    G2=G.copy()
    for i in G:
        dlist.append([i,G.degree(i)])
    d=sorted(dlist,key=itemgetter(1), reverse=True)
    #print(degreelist)
    newlist=d.copy()
    totalnbrs=[]
    for i in range(len(d)):                           #removing degree 0 vertices
        if d[i][1]==0:
            G2.remove_node(d[i][0])
            newlist.remove(d[i])
    d=newlist.copy()
    print(d)
    newlist=d.copy()
    print(S)
    if len(d)!=0:
        v=d[0][0]           #v is the highest degree vertex
        degv=d[0][1]
        if d[0][1]==2:
            for i in G2:
                for j in G2:
                    if i<j and G2.has_edge(i,j):
                        k=k-1
                        S.append(j)
        G=G2.copy()
        
        
        if v in G:
            Gd2=G.copy()
            vnbr=[n for n in Gd2.neighbors(v)]
            Gd2.remove_nodes_from(vnbr+[v])
            S1=S+vnbr
            return VertexCover(Gd2,k-degv,S1)
        if v in G:
            Gd1=G.copy()
            S2=S.append(v)
            Gd1.remove_node(v)
            
            return VertexCover(Gd1,k-1,S2)
          
    else:
        return True,S
 

        
"""

'''
def VertexCover(G,k,S):
    if G.edges()==[]:
        return True
    if k <=0:
       return False
    dlist=[]
    #G2=G.copy()
    for i in G:
        dlist.append([i,G.degree(i)])
    d=sorted(dlist,key=itemgetter(1), reverse=True)
    #print(degreelist)
    newlist=d.copy()
    totalnbrs=[]
    for i in range(len(d)):                           #removing degree 0 vertices
        if d[i][1]==0:
            G.remove_node(d[i][0])
            newlist.remove(d[i])
    d=newlist.copy()
    print(d)
    newlist=d.copy()
    #print(S)
    if len(d)!=0:
        v=d[0][0]           #v is the highest degree vertex
        degv=d[0][1]
        if d[0][1]<=2:
            for i in G:
                for j in G:
                    if i<j and G.has_edge(i,j):
                        return True

        Gd2=G.copy()
        Gd1=G.copy()
        vnbr=[n for n in Gd2.neighbors(v)]
        Gd2.remove_nodes_from(vnbr+[v])
        Gd1.remove_node(v)
        return VertexCover(Gd1,k-1,S+[v]) or VertexCover(Gd2,k-degv,S+vnbr)
    else:
        return True    
        



'''


def VertexCover(G,k,S):
    if G.edges()==[]:
        return S
    if k <=0:
       S=[]
       return S
    dlist=[]
    #G2=G.copy()
    for i in G:
        dlist.append([i,G.degree(i)])
    d=sorted(dlist,key=itemgetter(1), reverse=True)
    #print(degreelist)
    newlist=d.copy()
    totalnbrs=[]
    for i in range(len(d)):                           #removing degree 0 vertices
        if d[i][1]==0:
            G.remove_node(d[i][0])
            newlist.remove(d[i])
    d=newlist.copy()
    #print(d)
    newlist=d.copy()
    #print(S)
    if len(d)!=0:
        v=d[0][0]           #v is the highest degree vertex
        degv=d[0][1]
        if d[0][1]<=2:
            for i in G:
                for j in G:
                    if i<j and G.has_edge(i,j):
                        S=S+[j]
                        k=k-1
                        

        Gd2=G.copy()
        Gd1=G.copy()
        vnbr=[n for n in Gd2.neighbors(v)]
        Gd2.remove_nodes_from(vnbr+[v])
        Gd1.remove_node(v)
        S1=VertexCover(Gd1,k-1,S+[v]) 
        S2=VertexCover(Gd2,k-degv,S+vnbr)
        if len(S1)!=0 and len(S2)!=0:
            if len(S1)<len(S2):
                S=S1
            else:
                S=S2
        return S
    else:
        return S    
        

VC0, VI0=[],[]
ans=Generator(Gf, 0, VC0, VI0, nodes,[])
Slist=[]
print('Partitioning to cliques and independent sets:\n')
print('\n')
for i in ans:
    VC=i[0]
    VI=i[1]
    print(VC,VI)
    print('\n')
    GVC=nx.Graph()
    GVI=nx.Graph()
    for i in VC:
        for j in VC:
            if i !=j and Gf.has_edge(i,j):
                GVC.add_nodes_from([i,j])
                GVC.add_edge(i,j)
    for i in VI:
        for j in VI:
            if i !=j and Gf.has_edge(i,j):
                GVI.add_nodes_from([i,j])
                GVI.add_edge(i,j)
    Gnew=nx.Graph()
    GVCcom=nx.complement(GVC)
    Gnew.add_nodes_from(GVI.nodes())
    Gnew.add_nodes_from(GVC.nodes())
    for u in GVI:
        for v in GVI:
            if u!=v and GVI.has_edge(u,v):
                Gnew.add_edge(u,v)
    for u in GVCcom:
        for v in GVCcom:
            if u!=v and GVCcom.has_edge(u,v):
                Gnew.add_edge(u,v)  
    S=VertexCover(Gnew,k,Sstar)
    S1=[]
    for i in S:
        if i not in S1 and i in H.nodes():
            S1.append(i)
    Slist.append(S1)
    #print(Slist)
'''
for i in flaglist:
    if i<=k:
        print('YES Instance of Vertex Cover and Split vertex deletion.',i)
        sys.exit()
print('NO Instance')

'''
'''
for i in flaglist:
    if i==True:
        print('YES Instance of Vertex Cover and Split vertex deletion.',i)
        sys.exit()
print('NO Instance')
    
'''
print('Number of vertices in G=',H.number_of_nodes())
print('Number of edges in G=',H.number_of_edges())
print('k=',k)
print('Result:')
print('-------')
#Slist.sort(key=len)
for i in Slist:
    if len(i)<=k and len(i)!=0:
        print('YES Instance of Vertex Cover and Split vertex deletion.',i)
        sys.exit()
print('NO Instance')






















