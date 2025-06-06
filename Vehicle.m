function [] = Vehicle(center_location,L,H,theta,rgb)
center1=center_location(1);
center2=center_location(2);
R= ([cos(theta), -sin(theta); sin(theta), cos(theta)]);
X=([-L/2, L/2, L/2, -L/2]);
Y=([-H/2, -H/2, H/2, H/2]);
for i=1:4
T(:,i)=R*[X(i); Y(i)];
end
x_lower_left=center1+T(1,1);
x_lower_right=center1+T(1,2);
x_upper_right=center1+T(1,3);
x_upper_left=center1+T(1,4);
y_lower_left=center2+T(2,1);
y_lower_right=center2+T(2,2);
y_upper_right=center2+T(2,3);
y_upper_left=center2+T(2,4);
x_coor=[x_lower_left x_lower_right x_upper_right x_upper_left];
y_coor=[y_lower_left y_lower_right y_upper_right y_upper_left];
%patch('Vertices',[x_coor; y_coor]','Faces',[1 2 3 4],'Edgecolor',rgb,'Facecolor','none','Linewidth',1.2);
fill(x_coor, y_coor,rgb);
axis equal;
end