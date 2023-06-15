%%%%%Sismique réfraction %%%%%
%%Dromochroniques (x,y) %%%%%%
shots_file1=load('shots.dat');
sx=shots_file1(:,1); 
sz=shots_file1(:,2);
gx=shots_file1(:,3);
gz=shots_file1(:,4);
t_time=shots_file1(:,5);
offset=zeros(size(shots_file1,1),1);
clf
incr=0;

for i=1:24:size(shots_file1,1)
    x=(gx(i:i+23)-sx(i)+incr);y=t_time(i:i+23);
    x=[sx(i); x];y=[0; y];
    p = polyfit(x,y,6);
    f = polyval(p,x); 
    figure(1);plot(x,y,'x',x,f);hold on;
    incr=incr+20;
end 
hold off;
text(0,-0.003,'Sources');
ylabel('Time (s)');xlabel('Offset (x_{ij})');

%%%%Inverse généralisée %%%%%%%%%%%%%%%%%%%%%%%
%%Matrice L (mat_l) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mat_l=[nb_zeros_s,1,nb_zeros_gauche,1,nb_zeros_droite,Xij];
%nb_zeros_s :       Autant de zeros que de points de pas de tir
%nb_zeros_gauche :  Autant de zéros que de points d'offset
%nb_zeros_gauche :  Autant de zéros que de points entre le géophone...
%                   considéré et la fin du profil

%Xij :              Offset

s=1;g=1;
delta_s=0;delta_g=10;delta_tot=abs(shots_file1(1,1)-shots_file1(end,3))/delta_g;
m=size(shots_file1,1);
mat_l=zeros(m,delta_tot+2);
for i=1:m
    if (i>24 && mod(i,24)==1)
        delta_s=delta_s+2;
    end
    delta_p=abs(shots_file1(i,1)-shots_file1(i,3))/delta_g;
    nb_0_minus=delta_p-1;
    mat_nb_0_minus=zeros(1,nb_0_minus);
    nb_0_plus=(delta_tot-delta_s-nb_0_minus-s);
    mat_nb_0_plus=zeros(1,nb_0_plus);
    if delta_s==0
        mat_l(i,1:delta_tot+2)=[s,mat_nb_0_minus,g,mat_nb_0_plus,(delta_p*delta_g)];
    else
    mat_l(i,1:delta_tot+2)=[zeros(1,delta_s),s,mat_nb_0_minus,g,mat_nb_0_plus,(delta_p*delta_g)];
    end
end
figure(2);spy(mat_l);daspect([1 10 1]);
xlabel(['n+1 = ' num2str(size(mat_l,2))]);
ylabel(['m = ' num2str(size(mat_l,1))]);

%%Calcul de Vb %%%%%%%%%%%%%
mat_p=mldivide(mat_l,t_time);
mat_p_corr=mat_p(isnan(mat_p)==0);
vb=1/(mat_p(end));

%%Choix de Vw %%%%%%%%%%%%%
n_corr=size(mat_p_corr,1)-1;
disp(vb);
vw=[vb/3 vb/3+500 vb/3-500];
mat_z=zeros(n_corr,length(vw));
col=['r';'g';'b'];
% [C, ia, ic] = unique(gx,'sorted');
[C, ia, ic] = unique(gx);
elev=[sz(1);sz(25)/2;sz(25);sz(49)/2;sz(49);gz(ia)];
for j=1:length(vw)
    vw_temp=vw(j);
    mat_z(1:n_corr,j)=mat_p_corr(1:end-1)*(vb*vw_temp/(sqrt(vb^2-vw_temp^2)));
    figure(3);subplot(3,1,j);plot(elev-mat_z(1:n_corr,j));hold on;plot(elev);hold off;
    text(mat_z(j),1, (['v_{w}=' num2str(round(vw_temp)) ' m/s']), 'Color', col(j));
    if j==2
        ylabel('Elevation (m)');
    end
end

%%%%Statiques résiduelles %%%%%%%%%%%%%%%%
%%Correction des délais %%%%%%%%%%%%%%%%%%
for ik=1:3
    file=(['xc_data' num2str(ik) '.dat']);
    tr=load(file);
    amp_tr=reshape(tr(:,3),749,75);
    n_column=size(amp_tr,2);
    t_ref=sum(amp_tr(:,1:n_column),2);

    ik1=2*ik-1;
    figure(5);subplot(3,2,ik1);imagesc(tr(1:749,1),tr(1:749,2),amp_tr);
    title([file(1:2) file(4:end)])

    for ij=1:n_column
    
        [c_ww,lag] = xcorr(t_ref,amp_tr(:,ij),'coeff');
        [ind,I]=max(abs(c_ww));max_lag=lag(I);
        amp_tr(:,ij)=circshift(amp_tr(:,ij),max_lag);
    end
    
    ik2=ik*2;
    figure(5);subplot(3,2,ik2);imagesc(tr(1:749,1),tr(1:749,2),amp_tr);
    title('Après correction');
    
end   










