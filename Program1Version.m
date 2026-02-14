function Program1Version()

    clc; close all;

    epsilon1 = 2.5;
    mu1 = 0.01;

    epsilon2 = 5.0;
    mu2 = 0.01;

    E0_i = 2.0;
    omega = 2*pi;
    theta_deg = 45; 
    theta_i = deg2rad(theta_deg);

    L = 10;
    N = 300;
    z = linspace(-L, L, N);
    x = linspace(-L, L, N);
    [Z, X] = meshgrid(z, x);

    beta1 = omega * sqrt(epsilon1 * mu1);
    beta2 = omega * sqrt(epsilon2 * mu2);

    sin_theta_t = (beta1/beta2) * sin(theta_i);

    RT = false;
    if abs(sin_theta_t) > 1
        RT = true;
        theta_t = pi/2;
        disp('ATTENZIONE: Riflessione Totale');
    else
        theta_t = asin(sin_theta_t);
    end

    eta1 = sqrt(mu1/epsilon1);
    eta2 = sqrt(mu2/epsilon2);

    bx_i = beta1 * sin(theta_i);
    bz_i = beta1 * cos(theta_i);

    bx_r = bx_i;
    bz_r = -bz_i;

    bx_t = beta2 * sin(theta_t);
    bz_t = beta2 * cos(theta_t);

    if RT
        GammaTE = 1;
        TauTE = 0;
    else
        z1E = eta1/cos(theta_i);
        z2E = eta2/cos(theta_t);
        GammaTE = (z2E-z1E)/(z2E+z1E);
        TauTE = 1 + GammaTE;
    end

    E0_r = GammaTE * E0_i;
    E0_t = TauTE * E0_i;

    figure('Color', 'w', 'Name', 'Riflessione e Rifrazione TE', 'NumberTitle', 'off');
    colormap(redblue);
    
    ax = gca;
    hSurf = imagesc(z, x, zeros(size(Z))); % Placeholder per i dati
    hold on;
    
    xline(0, 'k-', 'LineWidth', 2);
    text(-L*0.9, L*0.9, 'Mezzo 1', 'FontSize', 12, 'FontWeight', 'bold');
    text(L*0.1, L*0.9, 'Mezzo 2', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w');
    
    axis xy; axis equal; axis tight;
    xlabel('Z'); ylabel('X');
    title(sprintf('Onda Piana: \\theta_i = %.1f°', theta_deg));
    colorbar;
    clim([-2 2]);

    t = 0;
    dt = 0.05;

    while ishandle(ax)
    
        E1 = E0_i * cos (omega * t - bx_i * X - bz_i * Z) + E0_r * cos (omega * t - bx_r * X - bz_r * Z);

        if RT
            az_t = sqrt((beta1^2) * (sin(theta_i)^2) - (beta2^2));
            decay = exp(-az_t * Z);
            E2 = E0_i * decay .* cos(omega * t - bx_t * X);
        else
            E2 = E0_t * cos(omega * t - bx_t * X - bz_t * Z);
        end

        E_tot = E1 .* (Z >= 0) + E2 .* (Z < 0);

        set(hSurf, 'CData', E_tot);

        t = t + dt;
        drawnow;
    end
end

function c = redblue()
    m = 100;
    r = (0:m-1)'/max(m-1,1); 
    c = [r, r, ones(m,1); ones(m,1), flipud(r), flipud(r)];
end
