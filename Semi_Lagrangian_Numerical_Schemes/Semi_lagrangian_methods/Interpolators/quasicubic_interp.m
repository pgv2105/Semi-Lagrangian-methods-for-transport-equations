

function z = quasicubic_interp(x,y,Xe,Ye,Ze,ind)

[n,m]=size(Xe);
%[ind,~] = knnsearch([Xe(:) Ye(:)],[x y],'K',4);
r = Xe(ind).^2+Ye(ind).^2;
ind=ind(find(r==min(r)));
[row,col] = ind2sub([n, m],ind);

% ind = [(row-1) (row-1); col (col+1)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z0 = interp1(Xe(st), Ze(st), x, 'linear');
% Y0 = Ye(st(1));

ind = [(row-1)*ones(1,4); (col-1):(col+2)];
st=sub2ind([n,m],ind(1,:),ind(2,:));
Z0 = interp1(Xe(st), Ze(st), x, 'linear');
Y0 = Ye(st(1));

ind = [(row)*ones(1,4); (col-1):(col+2)];
st=sub2ind([n,m],ind(1,:),ind(2,:));
Z1 = interp1(Xe(st), Ze(st), x, 'spline');
Y1 = Ye(st(1));

ind = [(row+1)*ones(1,4); (col-1):(col+2)];
st=sub2ind([n,m],ind(1,:),ind(2,:));
Z2 = interp1(Xe(st), Ze(st), x, 'spline');
Y2 = Ye(st(1));
 
ind = [(row+2)*ones(1,4); (col-1):(col+2)];
st=sub2ind([n,m],ind(1,:),ind(2,:));
Z3 = interp1(Xe(st), Ze(st), x, 'linear');
Y3 = Ye(st(1));
% 
% ind = [(row+2) (row+2); col (col+1)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z3 = interp1(Xe(st), Ze(st), x, 'linear');
% Y3 = Ye(st(1));

z = interp1([Y0 Y1 Y2 Y3], [Z0 Z1 Z2 Z3], y, 'spline');


end


% 
% function z = quasicubic_interp(x,y,Xe,Ye,Ze,ind)
% 
% [n,m]=size(Xe);
% %[ind,~] = knnsearch([Xe(:) Ye(:)],[x y],'K',4);
% r = Xe(ind).^2+Ye(ind).^2;
% ind=ind(find(r==min(r)));
% [row,col] = ind2sub([n, m],ind);
% 
% % ind = [(row-1) (row-1); col (col+1)];
% % st=sub2ind([n,m],ind(1,:),ind(2,:));
% % Z0 = interp1(Xe(st), Ze(st), x, 'linear');
% % Y0 = Ye(st(1));
% 
% ind = [(row-1)*ones(1,4); (col-1):(col+2)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z0 = interp1(Xe(st), Ze(st), x, 'linear');
% Y0 = Ye(st(1));
% 
% ind = [(row)*ones(1,4); (col-1):(col+2)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z1 = interp1(Xe(st), Ze(st), x, 'spline');
% Y1 = Ye(st(1));
% 
% ind = [(row+1)*ones(1,4); (col-1):(col+2)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z2 = interp1(Xe(st), Ze(st), x, 'spline');
% Y2 = Ye(st(1));
%  
% ind = [(row+2)*ones(1,4); (col-1):(col+2)];
% st=sub2ind([n,m],ind(1,:),ind(2,:));
% Z3 = interp1(Xe(st), Ze(st), x, 'linear');
% Y3 = Ye(st(1));
% 
% % ind = [(row+2) (row+2); col (col+1)];
% % st=sub2ind([n,m],ind(1,:),ind(2,:));
% % Z3 = interp1(Xe(st), Ze(st), x, 'linear');
% % Y3 = Ye(st(1));
% 
% z = interp1([Y0 Y1 Y2 Y3], [Z0 Z1 Z2 Z3], y, 'cubic');
% 
% 
% end