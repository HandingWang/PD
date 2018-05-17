% Authors:    Handing Wang, Yaochu Jin, Xin Yao
% University of Surrey, UK, and University of Birmingham, UK
% EMAIL:      wanghanding.patch@gmail.com, yaochu.jin@surrey.ac.uk, X.Yao@cs.bham.ac.uk
% WEBSITE:    http://www.surrey.ac.uk/cs/people/handing_wang/
% DATE:       March 2016
% ------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
% Handing Wang, Yaochu Jin, Xin Yao, Diversity Assessment in Many-Objective Optimization, Cybernetics, IEEE Transactions on, Accepted, 10.1109/TCYB.2016.2550502.
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
% ------------------------------------------------------------------------
function [ pd ] = PD( X )
% Usage: [ pd ] = PD( X )
%
% Input:
% X             -Objective values of the population n*m (n solutions with m objectives)
%
% Output: 
% pd            -PD value of population X
%
p=0.1;%lp norm setting
C=zeros(size(X,1),size(X,1));%connection array
D=zeros(size(X,1),size(X,1));%dissimilarity array
%Calculate the dissimilarity between each two solutions
for i=1:size(X,1)-1
    for j=i+1:size(X,1)
       d=sum(abs(X(j,:)-X(i,:)).^p,2).^(1/p);
       D(i,j)=d;
       D(j,i)=d;
    end
end
DMAX=max(max(D))+1;
D(logical(eye(size(D))))=DMAX;
n=size(X,1);
pd=0;

for k=1:n-1
    %Find the nearest neighbor to each solution according to D in each row.
    [d,J]=min(D,[],2);
    %Find solution i with the maximal di to its neighbor j
    [dmx,i]=max(d);

    while liantong(C,i,J(i))==1 %i and j are connected by previous assessed solutions
        if D(J(i),i)~=-1
            D(J(i),i)=DMAX; %Mark the connected subgraph
        end
        if D(i,J(i))~=-1
            D(i,J(i))=DMAX;
        end
        [d,J]=min(D,[],2);
        %Find solution i with the maximal di to its neighbor j
        [dmx,i]=max(d);
    end
    C(J(i),i)=1;
    C(i,J(i))=1;
    pd=pd+dmx;
    if D(J(i),i)~=-1
        D(J(i),i)=DMAX;%Mark the used dissimilarity di.
    end
    D(i,:)=-1;%Mark the chosen solution i
end
end

function [w]=liantong(C,I,J)
% Usage: [w]=liantong(C,I,J)
%
% Input:
% C             -Connection array
% I             -index I
% J             -index J
%
% Output: 
% w             -1 if solutions I and J are connected, 0 if solutions I and J are not connected.
%
V=I;
Child=find(C(V,:)==1);
if isempty(find(Child==J))==0
    w=1;
    return
else
    C(V,:)=0;
    C(:,V)=0;
    for i=1:size(Child,2)
        w=liantong(C,Child(i),J);
        if w==1
            return
        end
    end
end
w=0;
end


