% kldemonew.m
% Demonstration of Kernighan-Lin Bi-partitioning Algorithm
% (C) 2004 by Yu Hen Hu
% created: 9/4/2004
% updated: 9/11/2004 add a third example
% updated: 12/14/2004 add choice = 3, and interactive problem entry

clear all, close all,

% Circuit specification
% the circuit is specified as a n x n incident matrix c
% where n is even
% c(i,i) = 0, c(i,j) = c(j,i), 
% c(i,j) = 0 if there is no connection between 
%   vertex i to vertex j, and hence the cost of partition them is 0.
valid=0;
while ~valid,
   disp('Enter 0 to use circuit discussed in the class note (default);');
   disp('Enter 1 to use the circuit in Figure 2.2(b) of text book;');
   disp('Enter 2 to use another example;');
   disp('Enter 3 to use Figure 4.1(a) of textbook, p. 162;');
   disp('Enter 4 to enter connectivity matrix and initial paritition interactively:');
   chos=input('Enter a selection, ENTER key means default: ');
   if isempty(chos), chos=0; 
   end
   if ismember(chos,[0 1 2 3 4]), valid=1;
   else disp('Entry invalid, must be an integer 0, 1, 2, 3, or 4: ');
   end
end

if chos==0,
   c0=[0 1 0 0 0 0;
       1 0 1 1 0 0;
 	    0 1 0 0 0 0;
 	    0 1 0 0 1 1;
 	    0 0 0 1 0 1;
 	    0 0 0 1 1 0];
	   idxA0=[2 3 4]; idxB0=[1 5 6];
	   c0=c0([idxA0 idxB0],[idxA0 idxB0]); % permute c0 to fit initial partition
elseif chos==1,  
   c0=[0 1 2 3 2 4;
    	 1 0 1 4 2 1;
       2 1 0 3 2 1;
   	 3 4 3 0 4 3;
   	 2 2 2 4 0 2;
   	 4 1 1 3 2 0];
 		idxA0=[1 2 3]; idxB0=[4 5 6]; % initial partition node list in A and B
elseif chos==2, 
   c0=[0  1  0  4  2  0
       1  0  0  2  2  3
       0  0  0  0  6  0
       4  2  0  0  0  0
       2  2  6  0  0  4
       0  3  0  0  4  0];
    idxA0=[1 2 3]; idxB0=[4 5 6]; % initial partition node list in A and B
elseif chos==3,
   c0=[0 0 0 0 1 0 0 0 
       0 0 0 0 1 0 0 0
       0 0 0 0 0 1 0 0
       0 0 0 0 0 1 0 0
       1 1 0 0 0 0 1 0
       0 0 1 1 0 0 1 0
       0 0 0 0 1 1 0 1
       0 0 0 0 0 0 1 0];
   idxA0=[1 2 3 4]; idxB0 = [5 6 7 8];
elseif chos==4,
   c0=input('Enter c0 matrix, remember bracket: ');
   idxA0=input('Initial partition A (as a row vector) = ');
   idxB0=input('Initial partition B (as a row vector) = ');
else
   disp('Input invlid, re-enter between 0 to 3, or press return key: ');
end

% note that A and B are disjoint sets
disp('Initial partitions are:')
disp(['Partition A: ' int2str(idxA0)]);
disp(['Partition B: ' int2str(idxB0)]);
disp('Intial incident matrix [partition A; partition B] is:');
disp(c0);

c=c0; idxA=idxA0; idxB=idxB0; 
n=size(c0,1);   % c0 is a square, symmetric n x n matrix
n1=n/2; n2=n1;  % equal size bi-partition |A| = |B| = n/2

% a few words about indexing of nodes:
% the nodes are labeled as partition A: [1, 2, ..., n1=n/2]
%              partition B: [n1+1, n1+2,..., n]
% during iterations, the c matrix will reduce its size by 2
% since the two exchange nodes will be frozen during future iterations
% hence the index to the columns and rows of the c matrix will be
% different from those of the node indices.
% in this program, we will use node indices to indicate nodes
% and matrix indices to refer to matrix column and row numbers
% note that only at the first iteration, these two are identical.
% later, they are different!

iter=1; % iteration count
done=0; % Boolean condition on whether iteration is done

while ~done, % while not yet done,
   disp(['*** Iteration ' int2str(iter) ' ***']);
   if iter==1, % during the first iteration, compute initial D value
    % compute internal cost, external cost, D value and g value
    % this will become a callable m-file later
    % compute internal cost, 1 x n
    intcost=[sum(c(1:n1,1:n1)') sum(c(n1+1:n,n1+1:n)')]; 
    % compute external cost, 1 x n
    extcost=[sum(c(1:n1,n1+1:n)') sum(c(n1+1:n,1:n1)')]; 
    % compute D values, 1 x n
    Dvalue=extcost-intcost;  % first n1 is partition A, next n2 is partition B
    disp(['Initial partition cost = ' int2str(sum(extcost(1:n1)))]);
   end % otherwise, the D values will be updated later. 
   
   % compute g value to determine candidate exchange nodes
   % compute all combinations of partition A nodes and partition B
   %   nodes among nodes that are not fixed. 
   gmat = Dvalue(1:n1)'*ones(1,n2)+ones(n1,1)*Dvalue(n1+1:n)...
    -2*c(1:n1,n1+1:n);  % gmat is n1 x n2 that is of the same
   % dimension as c(1:n1,n1+1:n). The row indices of gmat is the 
   % same as those of the c matrix. The column indices of gmat 
   % are the same as those of the c matrix minus n1
   [mtmp,id1]=max(gmat);  % find max g of each column of gmat matrix
   % mtmp, id1 are both 1 x n2
   [g(iter),id2]=max(mtmp); % g(iter) is max g of gmat matrix
   % id2 is the column index of the max g of the gmat matrix
   ida=id1(id2); idb=n1+id2;  % c matrix indices of two exchange nodes
   % that yield maximum g (gain) in reducing cut set cost.
   disp('The g matrix is:')
   disp(gmat);
   
   % now convert into node indices for the two selected nodes
   bidx(iter) = idxB(id2); % node index for node b
   aidx(iter) = idxA(ida);  % node index for node a
   disp(['iter = ' int2str(iter) ', nodes to be exchanged: ' ...
     int2str(aidx(iter)) ', and ' int2str(bidx(iter)) ...
    ' Max. g = ' num2str(g(iter))]);
   % In summary, the maximum g value corrsponds to exchange 
   % node aidx(iter) in A and node bidx(iter) in B, and these nodes correspond to
   % the ida, and idb rows and columns in the current c matrix
   
   % Now consider if it is necessary to go to the next iteration
   % Here for demonstration purpose, we use the naive criteria that
   % when c matrix is reduced to a 2 x 2 matrix in this iteration, then
   % it is done. Otherwise, we proceed.
   
   if size(c,1)==2, done=1; % stop the next iteration, done
   else
    iter=iter+1; % move to the next iteration
    % First remove node aidx(iter) and bidx(iter) from the node indices of 
    % A and B, idxA and idxB respectively
    idxA=setdiff(idxA,aidx);  % A-{a}
    idxB=setdiff(idxB,bidx);  % B-{b}
    
    % next, remove the ida, idb rows and columns of the c matrix
    idv=[setdiff([1:n],[ida idb]) [ida idb]]; % permute the indices of c
    c1 = c(idv,idv); % move ida, idb rows and columns to the last two rows
    % and last two columns 
    % before removing the last two columns and rows, we may update the 
    % Dvalues of remain nodes in partitions A and B due to the selection of 
    % nodes aidx and bidx for exchange
    % First permute the Dvalue vector:
    Dvalue=Dvalue(idv(1:n-2)); % move the Dvalue of a and b to last two entries
    % and then drop them. Dvalue is now 1 x (n-2)
    % 
    % D' = D + 2*c(idx in the same partition) - 2*c(idx in different partition)
    % Hence, partitions A and B must be updated separately
    % (I) partition A 
    % D' = D + 2*c(a,vA) - 2*c(b,vA)
    % correspond to first n1-1 column/rows of the c1 matrix and
    % ida is the second to last row/column of the c1 matrix
    % hence c(a,vA) is c1(n-1,1:n1-1), c(b,vA) is c1(n, 1:n1-1)
    Dvalue(1:n1-1)=Dvalue(1:n1-1)+2*c1(n-1,1:n1-1)-2*c1(n,1:n1-1);
    % (II) partition B
    % correspond to n1:n-2 column/row of the c1 matrix
    % D' = D + 2*c(b,vB) - 2*c(a,vB)
    % c(a,vB) is c1(n-1,n1:n-2), c(b,VB) is c1(n,n1:n-2)
    Dvalue(n1:n-2)=Dvalue(n1:n-2)+2*c1(n,n1:n-2)-2*c1(n-1,n1:n-2);
    % Dvalue is reduced length by 2
    disp('Updated D values are:');
    disp(Dvalue)

    % Now, reduced c matrix with ida, idb rows and columns removed
    c= c1(1:n-2,1:n-2); 
    n = size(c,1); n1=n/2; n2=n1; % update c matrix indices
    % Note that idxA is now 1 x n1, idxB is 1 x n2;
    disp('Press any key to continue to the next iteration ...')
    pause
   end % if not terminate, update D'
end % of while loop

% Now that all the g values are computed, we need to determine k,
% the number of exchanges needed.
% k is the maximum partial sum of the g vector.
% triu(toeplitz(g)) gives an upper triangular matrix. The first column
% contains g(1), the second column has g(2), g(1), the 3rd column has
% g(3), g(2), g(1), etc.
% sum(triu(toeplitz(g))) gives a running partial sum of the g sequence.
% K is the index of the maximum entry of the running partial sum.
[tmp,K]=max(sum(triu(toeplitz(g))));
disp('Decision ...');
disp(['Exchange first ' int2str(K) ' pairs of nodes']);
disp('Final partition = ');
disp(['Partition A: ' int2str(union(setdiff(idxA0,aidx(1:K)),bidx(1:K)))]);
disp(['Partition B: ' int2str(union(setdiff(idxB0,bidx(1:K)),aidx(1:K)))]);

