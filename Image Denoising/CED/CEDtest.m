%Runs the Coherence Enhancing Diffusion filter for different discretization
%schemes

I = im2double(imread('sync_noise.png'));
JS = CoherenceFilter(I,struct('T',15,'rho',10,'Scheme','S'));
JN = CoherenceFilter(I,struct('T',15,'rho',10,'Scheme','N'));
JR = CoherenceFilter(I,struct('T',15,'rho',10,'Scheme','R'));
JI = CoherenceFilter(I,struct('T',15,'rho',10,'Scheme','I'));
JO = CoherenceFilter(I,struct('T',15,'rho',10,'Scheme','O'));
figure, 

subplot(2,3,1), imshow(I), title('Before Filtering');
subplot(2,3,2), imshow(JI), title('Standard Scheme');
subplot(2,3,3), imshow(JN), title('Non Negative Scheme');
subplot(2,3,4), imshow(JI), title('Implicit Scheme');
subplot(2,3,5), imshow(JR), title('Rotation Invariant Scheme');
subplot(2,3,6), imshow(JO), title('Optimized Scheme');