function [A1,A2,A3,grad_v,tangent_basis,aires,index1,index2,A2index1,A2index2,A3index1,A3index2]=...
    prepare_advection2(V,faces,vertices,tau,grad_v,tangent_basis,aires,index1,index2,A2index1,A2index2,A3index1,A3index2)
nc=size(vertices,1); % number of vertices
nt=size(faces,1); % number of faces
if nargin<5
    % first computation of geometrical quantities
    flag=1;
    [cont_grad_v,grad_v,tangent_basis,tg_basis,aires,norm_faces,norm_vertices,index1,index2]=matrices_contraintes_ter(faces,vertices,3);
    A2index1=0;A2index2=0;A3index1=0;A3index2=0;
else
    flag=0;
end
    
clear cont_grad_v tg_basis norm_faces norm_vertices

% Projection of grad(wi) on the tangent plane

P_grad_v_c=cell(1,3); 
for s=1:3
    for k=1:2
       P_grad_v_c{s}(:,k)=sum(grad_v{s}.*tangent_basis{k,s},2);
    end
    P_grad_v{s}=repmat(P_grad_v_c{s}(:,1),1,3).*tangent_basis{1,s}+repmat(P_grad_v_c{s}(:,2),1,3).*tangent_basis{2,s};
end

% Local matrices

A1=ui_uj(nc,faces,aires,index1,index2,tau);
[A2,A2index1,A2index2]=gu_i_v_uj(nc,nt,aires,faces,V,P_grad_v,flag,A2index1,A2index2);
[A3,A3index1,A3index2]=gui_v_guj_v(nc,nt,aires,faces,V,P_grad_v,flag,A3index1,A3index2,tau);

function [A1]=ui_uj(nc,faces,aires,index1,index2,tau) 
% matrix for the integrals SS u_iu_j dx
termes=(1/tau)*repmat(aires/12,1,3);  
sympart=(sparse(index1,index2,termes,nc,nc));
A1=sparse(faces,faces,repmat(2*(1/tau)*aires,1,3),nc,nc)+sympart+sympart'; 

function [A2,A2index1,A2index2]=gu_i_v_uj(nc,nt,aires,faces,V,P_grad_v,flag,A2index1,A2index2) 
% matrix for the integrals SS (grad(u_i).v)u_j dx
A2val=[];
if flag
A2index1=[];
A2index2=[];
end

for s1=1:3
    for s2=1:3
        A2val=[A2val,(aires/24).*dot(P_grad_v{s1},V(faces(:,1),:)+V(faces(:,2),:)+V(faces(:,3),:)+V(faces(:,s2),:),2)];
        if flag
        A2index1=[A2index1,faces(:,s1)];
        A2index2=[A2index2,faces(:,s2)];
        end
    end
end

A2=sparse(A2index1,A2index2,A2val,nc,nc);

function [A3,A3index1,A3index2]=gui_v_guj_v(nc,nt,aires,faces,V,P_grad_v,flag,A3index1,A3index2,tau) 
% matrice for the integrals SS (grad(u_i).v)(grad(u_j).v) dx
A3val=[];
if flag
A3index1=[];
A3index2=[];
end
for s1=1:3
    for s2=1:3
        for a=1:3
            for b=1:3
                A3val=[A3val,
                    tau*((a==b)/72+1/72)*aires.*dot(P_grad_v{s1},V(faces(:,a),:),2).*...
                    dot(P_grad_v{s2},V(faces(:,b),:),2)];
                if flag
                A3index1=[A3index1,faces(:,s1)];
                A3index2=[A3index2,faces(:,s2)];
                end
            end
        end
    end
end


A3=sparse(A3index1,A3index2,A3val,nc,nc);

