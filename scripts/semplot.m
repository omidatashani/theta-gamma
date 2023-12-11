
function semplot(input,t, win, c1, c2, c3)

for k = 1 : size(input,1)
    jnk(k,:) = ndass_smooth(input(k, :), win);
end
matr = jnk(:,:);
inte = 2.5;
y  = nanmean(matr);
y1 = nanmean(matr)+nanstd(matr)./(size(matr, 1).^.5);
y2 = nanmean(matr)-nanstd(matr)./(size(matr, 1).^.5);
y3 = y2(size(matr , 2):-1:1);
y2 = y3;
Y  = [y1  y2];

x1 = t;
x2 = sort(t , 'ascend');
X  = [x1  x2];
fill(X, Y, [c1/inte c2/inte c3/inte], 'LineStyle', 'none', FaceAlpha=0.25);
hold on;
plot(t, y,  'color', [c1 c2 c3],'LineWidth',1);