f_sampl= 100e6;

n_elem_x= 256;
dim_x= 1/1000;

dim_y= [2 2]/1000;
n_elem_y= length(dim_y);
w_x= dim_x/4;
w_y= 1/4000;
focus=[0 0 30]/1000;
e_aperture = xdc_linear_multirow (n_elem_x,dim_x, n_elem_y,dim_y,w_x,w_y,1,2, focus);
%% 
set_sampling(f_sampl);
data = xdc_get(e_aperture,'rect');
[n,m]=size(data);
%% 
colormap(hot(128));
for i=1:m
x=[data(11,i), data(20,i); data(14,i), data(17,i)]*1000;
y=[data(12,i), data(21,i); data(15,i), data(18,i)]*1000;
z=[data(13,i), data(22,i); data(16,i), data(19,i)]*1000;
c= 1;
hold on
surf(x,y,z,c)
end
% Put som axis legends on
Hc = colorbar;
view(3)
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
grid
axis('image')
hold off

