Nodes = 20;
Ktotal = 4; 
DV = Nodes*Nodes*(2+Ktotal) + Ktotal;
obj = zeros(1,DV);
for i=1:Nodes
    for j=1:Nodes
        for k=1:Ktotal
        obj(1,varindex_3(1,i,j,k,Nodes,Ktotal)) = 1;
        obj(1,varindex_3(2,i,j,k,Nodes,Ktotal)) = 2;
        obj(1,varindex_3(3,i,j,k,Nodes,Ktotal)) = 2+k;
        obj(1,varindex_3(4,i,j,k,Nodes,Ktotal)) = 7;
        end
    end
end

        
       