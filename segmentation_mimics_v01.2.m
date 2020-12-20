%CELL'S MIMICS SEGMENTATION CODE !!!!!!!! NO TRACKING !!!!!!! Will segment
%one image and apply the crop to the whole tiff stack. If your images are
%drifting, please apply an image registration pipeline first. 


%Define the location of your Tiff stack with GFP channel. Only one channel
%color stack is supportert right now.
Path = 'K:\Research\AG Niederholtmeyer\microscopy\20201216\single position export\Single channel GFP\C1-4-point-communication_pos1.tif'
Frame_to_analyse = 49; %Frame to analyse 
mimics_not_analyse = []
imageSizeX = 2044;
imageSizeY = 2048;


tiff_info = imfinfo(Path); % return tiff structure, one element per image
tiff_stack = imread(Path, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(Path, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

se = strel('disk',50); % Filter matrice. You may need to tune it according to your image size.
se2 = strel('disk',20); %Filter matrice. You may need to tune it according to your image size.

figure;

subplot(3,3,1)
I = tiff_stack(:,:,Frame_to_analyse);
I_eq2 = imadjust(I); % Contrast enhancement. Different methods will be available here.
imshow(I_eq2);  
title('Contrast enhanced')


subplot(3,3,2)
Io = imopen(I_eq2,se); % Erode dilatation step
imshow(Io)
title('Opening')

subplot(3,3,3)
Ie = imerode(I_eq2,se);
Iobr = imreconstruct(Ie,I_eq2);
imshow(Iobr)
title('Opening-by-Reconstruction')

subplot(3,3,4)
Iobrd = imdilate(Iobr,se2);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
imshow(Iobrcbr)
title('Reconstruction II')


I_eq = adapthisteq(Iobrcbr);
subplot(3,3,5)
imshow(I_eq);
title('Contrast enhanced II')

subplot(3,3,6)
bw = im2bw(I_eq, graythresh(I_eq));
imshow(bw)

bw2 = imfill(bw,'holes');
 

% se = ones(5,5)
se = strel('disk',20);

bw3 = imopen(bw2, se);



% Ie = imerode(bw3,se);
% Iobr = imreconstruct(Ie,bw3);

 
bw4 = bwareaopen(bw3, 1000);
subplot(3,3,7)
imshow(bw4)
% 
% bw5 = imerode(bw4, se);
%  imshow(bw5)
% 
% bw6 = imclearborder(bw4,4);
% subplot(2,3,6)
% imshow(bw6)

% 
% bw4_perim = bwperim(bw5);
% overlay1 = imoverlay(I_eq2, bw4_perim, [.3 1 .3]);
% subplot(2,3,6)
% 
% imshow(overlay1)

[centers,radii] = imfindcircles(I_eq2,[70 100],'Sensitivity',0.90,'Method','twostage')
subplot(3,3,8)
imshow(I_eq2)
h = viscircles(centers,radii-20);

%Display the labels of the found circles on the image
for k=1:length(radii)
         metricB1_string = sprintf('%2.0f',k);  
     text(centers(k,1)+60,centers(k,2)+60,metricB1_string,'color','w','HorizontalAlignment', 'center','VerticalAlignment', 'middle');
end


f = msgbox("Draw a circle in the background")

roi = drawcircle() %Ask the user for the background crop

yBack = roi.Center(2);
xBack = roi.Center(1);
RadiBack = roi.Radius;


%Analysis and plotting part 


Full_GFP=[] %contains the extracted fluorescent data
MeangreenBackground=[]



for k=1:size(tiff_stack,3)
I_i = tiff_stack(:,:,k);   %Update this one with the fluorescent image if you selected the bright-field for segmentation


[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
Meangreen = [];
n=1;

for i=1:length(radii)
    if ismember(i,mimics_not_analyse)
        continue; 
    end
centerX = centers(i,1);
centerY = centers(i,2);
radius = radii(i)-20;

circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;



indexMimics(n)=i;
Meangreen(n) = mean(I_i(circlePixels));

n=n+1;
end 
circlePixelsBackground = (rowsInImage -yBack).^2 ...
    + (columnsInImage - xBack).^2 <= RadiBack.^2;

MeangreenBackground(k) = mean(I_i(circlePixelsBackground));
Full_GFP(k,:)=Meangreen

end

 
figure,
hold all
colors = colormap(cbrewer('seq', 'Greens', 500));
colormap(colors);
for k=1:length(Full_GFP)
    plot(Full_GFP(:,k)'-MeangreenBackground,'Color',[0.4706 0.7765 0.4745])
end
    plot(mean(Full_GFP(:,:)'-MeangreenBackground),'k','LineWidth',2)
xlabel('Time (hours)')
ylabel('Fluorescence (A.U.)')
set(gca, 'XTick', 0:12:40);
set(gca, 'XTickLabel', {'0','1','2','3'});   
xlim([0 24])
ylim([0 inf])   


