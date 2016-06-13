# set wd
setwd("Z:/sandbox/cpProfile/")
matlab.script <- dir()[grep("\\.m", dir())]

# load function 
library(matconv)
cp_R <-mat2r(matlab.script)
str(cp_R)

# save to text file
file.create("cp_R2.R")
fileConn <- file("cp_R2.R")
writeLines(cp_R$rCode, fileConn)
close(fileConn)
file.show("cp_R2.R")

# create dummy data to test translated snipets
dummy.mat <- matrix(2,length(LETTERS),length(LETTERS),dimnames=list(LETTERS,LETTERS))
dummy.mat2 <- Matrix::Matrix(dummy.mat)
dummy.g <- igraph::graph.adjacency(dummy.mat)
dummy.net <- intergraph::asNetwork(dummy.g)

###
### START CODE
###

# matlab translated code
cpProfile <- function(net.object, directed=NULL) { } # leave function closed until translation is completed

# load required packages
require(Matrix)
require(network)
require(lattice)
require(igraph)

# create netname from matrix, igraph, or network objects and assign labels
if(class(net.object) %in% c("matrix", "dsyMatrix", "dscMatrix", "dsparseMatrix", "dsRMatrix", "dtCMatrix", "dtpMatrix", "dtRMatrix", "dtrMatrix")) { 
    print("Matrix object detected. Converting to sparse Matrix format")
    if(is.null(rownames(net.object))) { labels <- 1:nrow(net.object) }
    else { labels <- rownames(net.object) }
    netname <- Matrix::Matrix(net.object, dimnames=list(labels,labels)) }
if(class(net.object)=="igraph") {
    print("igraph object detected. Converting to sparse Matrix format")
    try(label.value <- grep("name", igraph::list.vertex.attributes(net.object), value = T)[1])
    if(length(label.value)==0 | is.null(label.value)) { labels <- 1:igraph::vcount(net.object) }
    else { labels <- igraph::get.vertex.attribute(net.object, label.value) }
    netname <- Matrix::Matrix(igraph::get.adjacency(net.object), dimnames=list(labels,labels)) }
if(class(net.object)=="network") {
    print("Network object detected. Converting to sparse Matrix format")
    try(label.value <- grep("name", network::list.vertex.attributes(net.object), value = T)[1])
    if(length(label.value)==0 | is.null(label.value)) { labels <- 1:network::network.size(net.object) }
    else { labels <- network::get.vertex.attribute(net.object, label.value) }
    netname <- Matrix::Matrix(igraph::get.adjacency(intergraph::asIgraph(net.object)), dimnames=list(labels,labels)) }

# check if convertion was sucesscul
if(!exists("netname")) { stop(print("cpProfile requires a matrix, network, or igraph object")) }

# start routine
if(any(netname>1)) { cat("Weighted network detected\n\n") }
cat("Profiling Core-Periphery\n\n")

# calculate algorithm metrics
A <- netname
k_in <- colSums(A)              # row vector of node in-weights (or in-degrees)
k_out <- rowSums(A)             # column vector of node out-weights (or out-degrees)
k_tot <- k_out+k_in             # total degree (or twice the degree, if undirected)
m <- sum(k_in)                  # total weight (or total number of links) in the network
N <- length(k_in)               # number of nodes
Abin <- as(A, "lgCMatrix")      # binary adjacency matrix
if(directed==T) {directed <- 1} # use info for directed network
if(directed==F) {directed <- 0} # use info for undirected network
if(is.null(directed)) {
    if(any(rowSums(A)==colSums(A))) {
        directed <- 0 }         # identify undirected network
    else {
        directed <- 1 }         # identify directed network
}
if(any(A>1)) { weighted <- 1 }  # identify weighted networks
if(!any(A>1)) { weighted <- 0}  # identify weighted networks

# save.image("cpProfile.Rdata")
# ## RESUME FROM HERE

# compute Markov matrix
cat("Computing the Markov matrix\n\n")
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


}




