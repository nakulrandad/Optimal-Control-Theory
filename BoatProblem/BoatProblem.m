%--------------------------------------------------------------------------
% BoatProblem.m
% The problem we try to solve here is the famous Zelmelo’s navigation
% problem. The goal is to minimize the time required by the ship to reach
% the desired position on the other bank from the start location.
%--------------------------------------------------------------------------
% Primary Contributor: Nakul Randad, Indian Institute of Technology Bombay
%--------------------------------------------------------------------------
clear all
global x0 y0 xf yf ty W V

% initial state
x0 = 0;
y0 = 0;

% end state
xf = 500;
yf = 150;

% Velocity of boat
V = 12; 

% Maximum velocity of riverwater
W = 4; 

% line distance for an initial time guess
tf0 = sqrt(xf^2+yf^2)/V;

% initial co-states
p0 = [-1; 1];

% solve the nonlinear equations
in = fsolve(@func,[p0;tf0]);

% start the animation
drawBoat(ty, W, [xf,yf]', in);

%%
function res = func(in)
    global x0 y0 xf yf ty norm_tf W V

    % rename variables
    p0 = in(1:2);
    tf = in(end);
    disp(in)
    norm_tf = tf; % Requried normalisation for numerical stability

    % time vector (transformation to fixed endtime)
    tspan = linspace(0, 1, 750);

    % integrate the canonical system
    [t, out] = ode45(@canon_sys,tspan,[x0;y0;p0]);

    % save and rename variables
    ty = [t out];
    x_tf = out(end,  1)' ;
    y_tf = out(end, 2)';
    p_tf = out(end,  3:4)' ;

    % evaluate the transversality condition 
    u_tf = atan(p_tf(2)/p_tf(1));
    H_tf = 1+p_tf(1)*V*cos(u_tf)+p_tf(2)*(4/(xf^2)*x_tf*(xf-x_tf)*W+V*sin(u_tf));

    % calculate the residuum
    res = [ x_tf-xf ;
            y_tf-yf;
           H_tf-0 ];
end


function outp = canon_sys(t,out)
    global norm_tf V W xf

    % rename variables
    p1 = out(3);
    p2 = out(4);
    x1 = out(1);

    % evaluate the stationary condition
    u = atan(p2/p1);

    % calculate the flow velocity
    d = ((4/xf)*x1-(4/xf^2)*x1^2)*W;

    % evaluate the canonical system
    outp = norm_tf*[ V*cos(u);...
                 d+V*sin(u);...
                 -W*p2*(4/xf-8/(xf^2)*x1);...
                 0 ];
end


function drawBoat(ty, W, sf, w)
    % position of the boat
    px = ty(:,2);
    py = ty(:,3);
    u = atan(ty(:,5)./ty(:,4));

    % create the river and drift vectors
    b = sf(1); % breadth 
    h = b/2; % height
    fig = figure('position', [50 500 600 500]);
    plot([0 0],[-h h],'k','LineWidth',4)
    title('TIME OPTIMAL RIVER CROSSING (Zermelos navigation problem)')
    hold on
    plot([b b],[-h h],'k','LineWidth',4)
    [x1,x2] = meshgrid(0:50:b,-h:100:h);
    q1 = 0*x1;
    q2 = 4/b^2*x1.*(b-x1)*W;
    hq = quiver(x1,x2,q1,q2);
    set (hq, "maxheadsize", 0.15);
    xlim([-100 b+100]);
    ylim([-h-100 h+150]);
    xlabel('x in m');
    ylabel('y in m');
    boatLength = 0.15*b;

    % plot trajectory start
    h0 = plot(0, 0, 'r');

    % plot start position
    h10 = plot([px(1)  px(1)+boatLength/6*cos(0)],[py(1) py(1)],'r','LineWidth',5);

    % plot boat at initial conditions
    h1 = plot([px(1)  px(1)+boatLength/2*cos(u(1))],[py(1) py(1)+boatLength/2*sin(u(1))],'k','LineWidth',5);

    % plot end position
    h2 = plot([b b+boatLength/6*cos(0)],[sf(2) sf(2)],'r','LineWidth',4);

    % plot some info
    h3 = text(-b*0.15, h*1.5, ['drift velocity: ', num2str(0)], 'fontsize', 11);
    h4 = text(-b*0.15, h*1.4, ['steering angle: ', num2str(0)], 'fontsize', 11);
    h5 = text(-b*0.15, h*1.3, ['trajectory time: ', num2str(0), ' s'], 'fontsize', 11);
    set(gca,'Color',[0.95 0.95 0.95],'XColor','k','YColor','k')
    set(gcf,'Color',[0.95 0.95 0.95])

    %start timer
    tic

    % animation in a for loop 
%     axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'animated_boat.gif';
    for i = 1:length(px)

        % draw trajectory
        set(h0, 'xdata', px(1:i),'ydata', py(1:i))

        % update pole and cart position
        set(h1, 'xdata', [px(i), px(i)+boatLength/2*cos(u(i))], 'ydata',[py(i), py(i)+boatLength/2*sin(u(i))]);

        % print drift velocity, angle and traj time
        set(h3, 'string',['drift velocity: ', num2str(4/b^2*px(i)*(b-px(i))*W), ' m/s']);
        set(h4, 'string',['steering angle: ', num2str(u(i)*180/pi), ' °']);
        set(h5, 'string',['trajectory time: ',num2str(w(3)*ty(i,1)),' s']);
        drawnow();

        %   % strop the animation when q is pressed
        %   if (kbhit (1) == 'q')
        %     break
        %   endif

        % pause the plotting
        pause(0.005);
        
%         % Capture the plot as an image 
%         frame = getframe(fig); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         if i == 1 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%         else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%         end 
        
    end

end