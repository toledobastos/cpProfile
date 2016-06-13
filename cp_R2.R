#
# !!! SIMPLIFIED VERSION FOR UNDIRECTED NETWORKS !!!
# C.P., Feb 2016
#
#######################################################################
#
#Profiling core-periphery structure in a (possibly) DIRECTED, WEIGHTED network.
#
#Please cite:
#F. Della Rossa, F. Dercole, C. Piccardi,
#Profiling core-periphery network structure by random walkers
#Scientific Reports, 3, 1467, 2013,
#http://dx.doi.org/10.1038/srep01467
#
#Copyright: 2013, Carlo Piccardi, Politecnico di Milano, Italy
#email carlo.piccardi@polimi.it
#
#Last updated: March 20, 2013
#
#######################################################################
#
#INPUT:
#The file
#   A_{netname}.mat
#must be in the working directory and must contain the following variables
#in Matlab binary format:
#   i) A : NxN weight matrix defining the (strongly connected) network.
#       A(i,j) is the weight of the link i->j.
#       If all the (nonzero) weights are 1, the network is actually UNWEIGHTED.
#       If A is symmetric, the network is actually UNDIRECTED.
#   ii) labels : (optional) 1xN cell vector of node labels (e.g., names)
#
#OUTPUT:
#On the screen:
#   i) the plot of the core-periphery profile;
#   ii) the cp-centralization C;
#   iii) the list of node labels with the corresponding coreness alpha_k
#        (ranked by the value of alpha_k)

clear all
close all
set(0,'Units','pixels')
scn <- get(0,'ScreenSize')

#######################################################################

#####name of the network: the file A_{netname}.mat will be loaded
#####UNCOMMENT the name of the network to be loaded

#-----the following networks are analyzed in the paper
#-----(see reference above): the datafiles are available in the
#-----distribution package

netname <- 'ff_und_bin'

#################################################################

#LOADING DATA, AND COMPUTING BASIC STATISTICS

disp([' '])
disp(['PROFILING CORE-PERIPHERY'])

#loads the NxN network matrix A
#and (optionally) a Nx1 cell "labels" containing label strings
load(strcat('A_',netname,'.mat'))
# A=full(A);
A <- sparse(A)
# ###A=double(A>0); #!!!!!!!!!!!!to test the binary net

#if labels do not exist in the uploaded file,
#creates fictitious labels which are simply the node numbers
if (length(find(char(who('*'))=='b'))==0){#labels do not exists in the file uploaded
    labels <- cell(length(A),1)
    for (i in 1:length(A)){
        labels(i) <- cellstr(num2str(i))
    }
}

#OPTIONAL: In the core-periphery profile algorithm, when many nodes
#attain the min, the one with smallest index is taken. If you want
#instead to randomize the selection, you can shuffle the node numbering by
#uncommenting the following block.

# N=length(A);
# rp=randperm(N);
# labels=labels(rp);
# RP=zeros(N);
# for i=1:N
#     RP(i,rp(i))=1;
# end;
# A=RP*A*RP^(-1);

disp(['Network: ',netname,' - N <- ',int2str(length(A))])

k_in <- sum(A)#row vector of node in-weights (or in-degrees)
k_out <- sum(A')'#column vector of node out-weights (or out-degrees)
k_tot <- k_out+k_in'#total degree (or twice the degree, if undirected)
m <- sum(k_in)#total weight (or total number of links) in the network
N <- length(k_in)#number of nodes

Abin <- double(A>0)#binary adjacency matrix
directed=sum(sum(Abin==Abin'))<N^2#directed=1 for directed networks
weighted=sum(sum(A==Abin))<N^2#weighted=1 for weighted networks

# # # disp(['Computing the Markov matrix...'])
# # # #creating the Markov matrix by row-normalizing A
# # # P=zeros(N,N);
# # # P=(diag(1./k_out))*A;
# # #
# # # disp(['Computing Markov chain asymptotic distribution...'])
# # # #computing Markov asymptotic distribution (x)
# # # if directed
# # #     AAA=eye(N)-P';
# # #     AAA(N,:)=1;
# # #     bbb=zeros(N,1);
# # #     bbb(N)=1;
# # #     x=AAA\bbb;
# # # else
# # #     x=k_in/sum(k_in);
# # # end;
# # #
# # # xP=diag(x)*P; #[xP]_ij = x_i * p_ij

# A=full(A);

#################################################################

#sorting nodes according to total degree
[L,nodelist] <- sort(k_in+k_out')
#current periphery and core:
#starting periphery from the least connected node
periph <- [nodelist(1)]
core <- setdiff([1:N],periph)){
#alpha_tmp(i) is the persistence probability of the periphery
#after the i-th node has been added
alpha_tmp <- zeros(1,N)

# # # #at each cycle, introducing the node that yields the
# # # #smallest increase in the pers.prob. of the periphery
# # # x_sum=sum(x(periph));
# # # xP_sum=sum(sum(xP(periph,periph)));
s_sum <- k_out(periph)
w_sum <- 0

tic#tracing the CPU time

for (i in 2:N-1){
    if (rem(i,100)==0){
        disp(['...adding node ',int2str(i),' of ',int2str(N)])
    }
    utest <- zeros(1,length(core))
    for (j in 1:length(core)){
        #computing the pers.prob. if node j is adedd to periphery

        utest(j) <- (w_sum+2*sum(A(periph,core(j)))+A(core(j),core(j)))/...
                 (s_sum+k_out(core(j)))
# # #         utest(j)=(xP_sum+sum(xP(core(j),periph))+sum(xP(periph,core(j)))+...
# # #                   xP(core(j),core(j)))/(x_sum+x(core(j)));
    }
    [uuu,jjj] <- sort(utest)
    alpha_tmp(i) <- uuu(1)
    #among the core nodes yielding minimal increase in the pers.prob.,
    #select the one with smallest total degree
    listmin=core(jjj(uuu==min(uuu)))
    k_listmin <- k_tot(listmin)'
    [kkk,lll] <- min(k_listmin)
    newnode <- listmin(lll)

    s_sum <- s_sum+k_out(newnode)
    w_sum <- w_sum+2*sum(A(periph,newnode))+A(newnode, newnode)

    periph <- [periph newnode]
    core <- setdiff([1:N],periph)){
# # #     x_sum=sum(x(periph));
# # #     xP_sum=sum(sum(xP(periph,periph)));

}

ttt <- toc#tracing the CPU time

#final step: the current periphery eventually includes the whole network
alpha_tmp(N) <- 1
periph <- [periph core]

#plotting the core-periphery profile
figure('OuterPosition',[1 1*scn(4)/2 1*scn(3)/3 scn(4)/2])
f1 <- get(0,'CurrentFigure')
figure(f1)
plot(1:N,alpha_tmp,'rx-',1:N,(0:N-1)/(N-1))
ylabel('core-periphery profile \it{\alpha_{k}}')
xlabel('number of nodes of \it{P_k}')
grid on

C <- sum((0:N-1)/(N-1)-alpha_tmp)/((N-2)/2)
disp([' '])
disp(['cp-centralization C <- ',num2str(C)])

#for export purposes: alpha(i) is the coreness of node i
for (i in 1:N){
    alpha(i)=alpha_tmp(periph==i)
}

disp(['CPU time (main cycle only) [sec] <- ',num2str(ttt)])
disp([' '])
disp(['Press a key to display node coreness...'])
pause
disp(['rank     node label     id     coreness alpha_k'])
for (i in N:-1:1){
#     disp([int2str(N-i+1),'   ',char(labels(periph(i))),'   ',char(idc(periph(i))),'   ',num2str(alpha_tmp(i))])
    disp([int2str(N-i+1),'   ',char(labels(periph(i))),'   ',num2str(alpha_tmp(i))])
}
disp(['rank     node label     coreness alpha_k'])

#saving the whole workspace
save(strcat('WksCP_',netname,'.mat'))

figure
plot(k_in,alpha,'o')
