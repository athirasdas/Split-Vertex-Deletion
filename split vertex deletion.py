import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
from scipy.io import mmread
import networkx.algorithms.isomorphism as iso
import math

import pandas as pd
Matrix = (mmread('input3.mtx'))
B = Matrix.todense()
#print(B)
df = pd.DataFrame(B, range(1, B.shape[0] + 1), range(1, B.shape[1] + 1))
#print(df)
G = nx.from_pandas_adjacency(df)
print(G)
#print(G.size())
#print(G.nodes())
print(G.nodes())

k=100
H=G.copy()
D=[]



class SVD:
    def __init__(self,G,k,D):
        self.G=G
        self.k=k
        self.D=D
    def RR1(self):
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
        #RR1 : to remove forbidden structures                                                                 
                                                                                
                                                                                #Cycle-5 is found and delete the nodes from G

        GM2 = iso.GraphMatcher(self.G, cycle5_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM2.subgraph_is_isomorphic():
            print('Cycle with 5 vertices found')
            dic2=GM2.mapping
            #print(dic1)
            D2.append(list(dic2.keys()))
            #print(D2)
            self.G.remove_nodes_from(D2[j])
            j=j+1
            #print(G.nodes())


                                                                                #Cycle-4 is found and delete the nodes from G
        GM1 = iso.GraphMatcher(self.G, cycle4_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM1.subgraph_is_isomorphic():
            print('Cycle with 4 vertices found')
            dic1=GM1.mapping
            #print(dic1)
            D1.append(list(dic1.keys()))
            #print(D1)
            self.G.remove_nodes_from(D1[i])
            i=i+1
            #print(G.nodes())


                                                                                #2K2 is found and delete the nodes from G    
        GM3 = iso.GraphMatcher(self.G, twok2_graph, node_match=iso.categorical_node_match(' ',' '), edge_match=iso.categorical_edge_match('', ''))
        while GM3.subgraph_is_isomorphic():
            print('2K2 is found')
            dic3=GM3.mapping
            #print(dic3)
            D3.append(list(dic3.keys()))
            #print(D3)
            self.G.remove_nodes_from(D3[l])
            l=l+1

        D4=D1+D2+D3
        self.D = self.D + [item for sublist in D4 for item in sublist]
        self.D.sort()
        print('D contain nodes:',self.D)    
        print('G-D contain nodes:',self.G.nodes())
        print('RR1 completed!!!!!!!!!')                                                                        #if |D|>5k it is a NO instance
        if (len(D)>(5*self.k)):
            print('NO Instance')
            return 
        else:
            return (self.G,self.k,self.D)
    
    def RR2(self):

        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        #print(cliques)


        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                Cstar=i
        print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        print('I*=',Istar)                                                       #Istar is I* ===> C* \disjoint_union\ I* is a split graph

        #RR2 begins

        V=self.G.nodes()
        E=self.G.edges()
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
            self.G.remove_nodes_from(vinter[0])
            self.k=self.k-1
            print(self.G.nodes)    

        print("RR2 completed")
        return (self.G,self.k,self.D)
        
    
    def RR3(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                Cstar=i
        print('C*=',Cstar)                                                                 # Cstar is C*
        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        print('I*=',Istar)                                                       #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=self.G.nodes()
        E=self.G.edges()
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
        print('RR3 completed!!!!')
        return (self.G,self.k,self.D)


    def RR4(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                Cstar=i
        print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        print('I*=',Istar)                                                       #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=self.G.nodes()
        E=self.G.edges()
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
        print('RR4 completed!!!')
        return (self.G,self.k,self.D)


    def RR5(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                Cstar=i
        print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        print('I*=',Istar)                                                       #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=self.G.nodes()
        E=self.G.edges()
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
        print('RR5 Completed!!')
        return (self.G,self.k,self.D)



    def RR6(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                Cstar=i
        print('C*=',Cstar)                                                                 # Cstar is C*

        Istar=[n for n in self.G.nodes() if n not in Cstar]                                     
        print('I*=',Istar)                                                       #Istar is I* ===> C* \disjoint_union\ I* is a split graph
        V=self.G.nodes()
        E=self.G.edges()
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
        
        print('RR6 Completed!!')
        return (self.G,self.k,self.D)


# #When none of the reduction rules apply, |C1*| ≤ 2k + 2 or |I1*| ≤ 2k + 2.

# #When none of the reduction rules apply, the number of vertices in C1* ∪ I1* is O(k^2).

# #When none of the reduction rules apply, the sets C* and I* contain O(k^3) vertices.



    def RR7(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                C=i
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
            for v in SY:
                self.G.remove_node(v)
                self.k=self.k-1
        print('RR7 completed!!!',S)
        return (self.G,self.k,self.D)



    def RR8(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                C=i
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
        print('RR8 completed!!!')
        return (self.G,self.k,self.D)




    def RR9(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                C=i
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
        print('RR9 completed!!!')
        return (self.G,self.k,self.D)



    def RR10(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                C=i
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
        print('RR10 completed!!!')
        return (self.G,self.k,self.D)



    def RR11(self):
        cliques=list(nx.find_cliques(self.G))                                          #returns the maximal cliques
        maxcliquelen=0                                                            #find a clique with maximum clique number 
        for i in cliques:
            if len(i)>maxcliquelen:
                maxcliquelen=len(i)
                C=i
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

        IYIZ=[n for n in IY if n in IZ]
        if len(IYIZ) > k+2:
            for v in IYIZ:
                for u in CfixUSC:
                    G.add_edge(v,u)
            for i in IYIZ:
                if len(IYIZ)==k+2:
                    break
                IYIZ.remove(i)
        print('RR11 completed!!!')
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


print(G11.nodes())
















































