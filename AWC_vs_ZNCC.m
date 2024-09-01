clear all
clc

% Skripta za optimizovano poređenje AWC i ZNCC modela sa unapređenim performansama

% Postavi random seme za reprodukcibilnost
rng(0);

% Dimenzije slike
imageSize = [512, 512];

% Generisanje osnovne slike (referentna slika)
I1 = imnoise(ones(imageSize), 'speckle', 0.1);

% Parametri za talasnu transformaciju
waveletName = 'db4';  % Diskretna talasna funkcija (Daubechies 4)
numLevels = 3;  % Broj dekompozicionih nivoa

% Parametri za kombinovanu korelaciju
alpha = 0.6;  % Težinski faktor za faznu korelaciju
beta = 0.4;   % Težinski faktor za amplitudnu korelaciju

% Definisanje nivoa šuma i deformacije
noiseLevels = linspace(0.01, 20, 5); % Nivoi Lorencovog šuma
deformationLevels = linspace(1, 10, 5); % Intenzitet deformacija

% Inicijalizacija promenljivih za čuvanje rezultata
awcmResults = zeros(length(noiseLevels), length(deformationLevels)); % Čuvanje stepena korelacije za AWC
znccResults = zeros(length(noiseLevels), length(deformationLevels)); % Čuvanje stepena korelacije za ZNCC
executionTimesAWC = zeros(length(noiseLevels), length(deformationLevels));
executionTimesZNCC = zeros(length(noiseLevels), length(deformationLevels));
ssimAWC = zeros(length(noiseLevels), length(deformationLevels));
ssimZNCC = zeros(length(noiseLevels), length(deformationLevels));
mseAWC = zeros(length(noiseLevels), length(deformationLevels));
mseZNCC = zeros(length(noiseLevels), length(deformationLevels));
pearsonAWC = zeros(length(noiseLevels), length(deformationLevels));
pearsonZNCC = zeros(length(noiseLevels), length(deformationLevels));

% Priprema za čuvanje slika za prikaz
savedImages = cell(length(noiseLevels), length(deformationLevels));

% Petlja kroz različite nivoe šuma i deformacija
for i = 1:length(noiseLevels)
    noiseLevel = noiseLevels(i);
    
    for j = 1:length(deformationLevels)
        deformationLevel = deformationLevels(j);

        % Generisanje deformisane slike sa Lorencovim šumom
        I2 = circshift(I1, round([5, 5] * deformationLevel)); % Zaokruživanje na ceo broj
        I2 = I2 + lorentzNoise(imageSize, noiseLevel);
        I2 = imrotate(I2, deformationLevel * 5, 'bilinear', 'crop'); % Rotacija
        I2 = imresize(I2, 1 + deformationLevel * 0.05); % Skaliranje
        I2 = imresize(I2, imageSize); % Vraćanje na originalne dimenzije

        % Čuvanje slike za prikaz
        savedImages{i, j} = I2;

        % Merenje vremena izvršavanja optimizovanog AWC modela
        tic;
        [cfs1, ~] = wavedec2(I1, numLevels, waveletName);
        [cfs2, ~] = wavedec2(I2, numLevels, waveletName);

        % Fazna korelacija
        phaseCorr = abs(sum(conj(cfs1) .* cfs2)) ./ (sqrt(sum(abs(cfs1).^2)) .* sqrt(sum(abs(cfs2).^2)));

        % Amplitudna korelacija
        amplitudeCorr = abs(cfs1) .* abs(cfs2);

        % Kombinovana korelacija
        combinedCorr = alpha * phaseCorr + beta * amplitudeCorr;

        % Normalizacija
        combinedCorr = combinedCorr / numLevels;
        executionTimeAWC = toc;
        executionTimesAWC(i, j) = executionTimeAWC;

        % Čuvanje rezultata za AWC
        awcmResults(i, j) = max(combinedCorr(:));

        % Merenje vremena izvršavanja ZNCC modela
        tic;
        zncc = normxcorr2(I1, I2);
        executionTimeZNCC = toc;
        executionTimesZNCC(i, j) = executionTimeZNCC;

        % Čuvanje rezultata za ZNCC
        znccResults(i, j) = max(zncc(:));

        % Izračunavanje SSIM vrednosti
        [ssimValAWC, ~] = ssim(imresize(combinedCorr, size(I1)), I1);
        [ssimValZNCC, ~] = ssim(imresize(zncc, size(I1)), I1);
        ssimAWC(i, j) = ssimValAWC;
        ssimZNCC(i, j) = ssimValZNCC;
        
        % Izračunavanje MSE vrednosti
        mseValAWC = immse(imresize(combinedCorr, size(I1)), I1);
        mseValZNCC = immse(imresize(zncc, size(I1)), I1);
        mseAWC(i, j) = mseValAWC;
        mseZNCC(i, j) = mseValZNCC;
        
        % Izračunavanje Pearsonovog koeficijenta korelacije
        pearsonValAWC = corr2(imresize(combinedCorr, size(I1)), I1);
        pearsonValZNCC = corr2(imresize(zncc, size(I1)), I1);
        pearsonAWC(i, j) = pearsonValAWC;
        pearsonZNCC(i, j) = pearsonValZNCC;

        % Prikaz rezultata za trenutni nivo šuma i deformacije
        fprintf('Noise Level: %.2f, Deformation Level: %.2f\n', noiseLevel, deformationLevel);
        fprintf('AWC: SSIM=%.4f, MSE=%.4f, Pearson=%.4f\n', ssimValAWC, mseValAWC, pearsonValAWC);
        fprintf('ZNCC: SSIM=%.4f, MSE=%.4f, Pearson=%.4f\n', ssimValZNCC, mseValZNCC, pearsonValZNCC);
        fprintf('Execution Time AWC: %.4fs, ZNCC: %.4fs\n\n', executionTimeAWC, executionTimeZNCC);
    end
end

% Prikaz grafikona poboljšanja
figure;
plot(noiseLevels, mean(awcmResults - znccResults, 2), '-o');
xlabel('Nivo Lorencovog šuma');
ylabel('Poboljšanje u %');
title('Poboljšanje optimizovanog AWC modela u odnosu na ZNCC');
grid on;

% Prikaz grafikona vremena izvršavanja
figure;
plot(noiseLevels, mean(executionTimesAWC, 2), '-o', noiseLevels, mean(executionTimesZNCC, 2), '-x');
xlabel('Nivo Lorencovog šuma');
ylabel('Vreme izvršavanja (s)');
legend({'AWC', 'ZNCC'}, 'Location', 'Best');
title('Vreme izvršavanja optimizovanog AWC vs ZNCC');
grid on;

% Prikaz uporednog SSIM grafikona
figure;
plot(noiseLevels, mean(ssimAWC, 2), '-o', noiseLevels, mean(ssimZNCC, 2), '-x');
xlabel('Nivo Lorencovog šuma');
ylabel('SSIM');
legend({'AWC', 'ZNCC'}, 'Location', 'Best');
title('Uporedni SSIM za optimizovani AWC i ZNCC model');
grid on;

% Prikaz uporednog MSE grafikona
figure;
plot(noiseLevels, mean(mseAWC, 2), '-o', noiseLevels, mean(mseZNCC, 2), '-x');
xlabel('Nivo Lorencovog šuma');
ylabel('MSE');
legend({'AWC', 'ZNCC'}, 'Location', 'Best');
title('Uporedni MSE za optimizovani AWC i ZNCC model');
grid on;

% Prikaz uporednog Pearson grafikona
figure;
plot(noiseLevels, mean(pearsonAWC, 2), '-o', noiseLevels, mean(pearsonZNCC, 2), '-x');
xlabel('Nivo Lorencovog šuma');
ylabel('Pearsonov koeficijent');
legend({'AWC', 'ZNCC'}, 'Location', 'Best');
title('Uporedni Pearsonov koeficijent za optimizovani AWC i ZNCC model');
grid on;

% Prikaz linijskih grafikona za AWC model u zavisnosti od Deformation Level i Noise Level
figure;
for j = 1:length(deformationLevels)
    plot(noiseLevels, awcmResults(:, j), '-o', 'DisplayName', ['Deformation Level ' num2str(deformationLevels(j))]);
    hold on;
end
xlabel('Noise Level');
ylabel('AWC Correlation');
legend('show');
title('AWC model: Korelacija u zavisnosti od Deformation Level i Noise Level');
grid on;

% Prikaz linijskih grafikona za ZNCC model u zavisnosti od Deformation Level i Noise Level
figure;
for j = 1:length(deformationLevels)
    plot(noiseLevels, znccResults(:, j), '-o', 'DisplayName', ['Deformation Level ' num2str(deformationLevels(j))]);
    hold on;
end
xlabel('Noise Level');
ylabel('ZNCC Correlation');
legend('show');
title('ZNCC model: Korelacija u zavisnosti od Deformation Level i Noise Level');
grid on;

% Prikaz devet karakterističnih slika korišćenjem IMSHOW
for i = 1:3
    for j = 1:3
        figure;
        imshow(savedImages{i, j}, []);
        title(sprintf('Noise: %.2f, Deformation: %.2f', noiseLevels(i), deformationLevels(j)));
    end
end

% Funkcija za generisanje Lorencovog šuma
function noise = lorentzNoise(imageSize, noiseLevel)
    % Parametri Lorencovog sistema
    sigma = 10;
    rho = 28;
    beta = 8/3;
    
    % Inicijalizacija
    dt = 0.01;
    numSteps = prod(imageSize);
    x = zeros(1, numSteps);
    y = zeros(1, numSteps);
    z = zeros(1, numSteps);
    
    % Početni uslovi (nasumični)
    x(1) = rand() * noiseLevel;
    y(1) = rand() * noiseLevel;
    z(1) = rand() * noiseLevel;
    
    % Integracija Lorencovog sistema
    for i = 2:numSteps
        dx = sigma * (y(i-1) - x(i-1)) * dt;
        dy = (x(i-1) * (rho - z(i-1)) - y(i-1)) * dt;
        dz = (x(i-1) * y(i-1) - beta * z(i-1)) * dt;
        x(i) = x(i-1) + dx;
        y(i) = y(i-1) + dy;
        z(i) = z(i-1) + dz;
    end
    
    % Normalizacija i generisanje slike šuma
    noise = reshape(z, imageSize);
    noise = noise / max(abs(noise(:))) * noiseLevel;  % Normalizuj šum na određeni nivo
end
