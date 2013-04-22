clear
fid=fopen('1.txt','r');
A=fscanf(fid,'%f',[7,inf]);
r(:,:)=A(2:4,:);
num=size(r(1,:),2);
dist(1:num)=999999999.0;
mindist(1:num)=0.0;
for i=1:num
	dist(i)=9999999.0;
	for j=1:num
		if ( j~=i )
			dist(j)=(r(1,i)-r(1,j)).^2+(r(2,i)-r(2,j)).^2+(r(3,i)-r(3,j)).^2;
		end
	end
	mindist(i) = min(dist);
end
[ss,idx]=min(mindist)
