% display_optical.m

image_disp = duw;
image_disp = repmat((image_disp-gmin)/(gmax-gmin),[1,1,3]);
image_disp = image_disp.^(1 - beta); % Brighten the image
imshow(real(image_disp),'XData',xp/1e3,'YData',yp/1e3);

axis xy
axis on
