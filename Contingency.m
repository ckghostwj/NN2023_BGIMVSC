function Cont=Contingency(Mem1,Mem2)
%CONTINGENCY Form contigency matrix for two vectors
% C=Contingency(Mem1,Mem2) returns contingency matrix for two
% column vectors Mem1, Mem2. These define which cluster each entity 
% has been assigned to.
%
% See also RANDINDEX.
%

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

%原来下一行代码中是不需要+1的，但handwritten中样本标签是0-9，其他数据集是从1开始。
%不+1的话下一行会出错。这个改动，只针对handwritten有效，其他数据集不用+1
Cont=zeros(max(Mem1) + 1,max(Mem2));

for i = 1:length(Mem1);
    try
        Cont(Mem1(i)+1,Mem2(i))=Cont(Mem1(i)+1,Mem2(i))+1;
    catch
        a=1;
    end
end
