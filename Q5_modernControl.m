A = [-1,0,-2,-2 ;0,-1,2,2 ;0,0,-1,-0.5 ;0,0,0,1]; % State matrix
B = [0;1;0;0]; % Control input matrix
C = [0, 0, 2, 1]; % Output matrix
D = zeros(2, 1);
n = size(A, 1);
obs_mat = [C; C*A];
for i = 2:n-1
    obs_mat = [obs_mat; C*A^i];
end
observability = rank(obs_mat) == n;

cont_mat = B;
for i = 2:n-1
    cont_mat = [cont_mat A^(i-1)*B];
end
controllability = rank(cont_mat) == n;

disp('Observability:');
disp(observability);
disp('Controllability:');
disp(controllability);

% Calculating eigenvalues of A
eigvals_A = eig(A);

% Checking unobservable modes
n = size(A, 1);
obs_mat = [C; C*A];
for i = 2:n-1
    obs_mat = [obs_mat; C*A^i];
end
unobservable_modes = [];
for i = 1:length(eigvals_A)
    if rank([eigvals_A(i)*eye(n) - A; obs_mat]) ~= n
        unobservable_modes = [unobservable_modes eigvals_A(i)];
    end
end

% Checking uncontrollable modes
cont_mat = B;
for i = 2:n-1
    cont_mat = [cont_mat A^(i-1)*B];
end
uncontrollable_modes = [];
for i = 1:length(eigvals_A)
    if rank([eigvals_A(i)*eye(n) - A cont_mat]) ~= n
        uncontrollable_modes = [uncontrollable_modes eigvals_A(i)];
    end
end

% results
disp('Unobservable Modes:');
disp(unobservable_modes);
disp('Uncontrollable Modes:');
disp(uncontrollable_modes);


% Calculate controllable and observable subspaces
n = size(A, 1);
cont_mat = B;
obs_mat = [C; C*A];
for i = 2:n-1
    cont_mat = [cont_mat A^(i-1)*B];
    obs_mat = [obs_mat; C*A^i];
end

% Perform QR decomposition to obtain controllable and observable subspaces
[Qc, ~] = qr(cont_mat);
[Qo, ~] = qr(obs_mat');

% Calculate uncontrollable and unobservable subspaces
uncontrollable_subspace = Qc(:, rank(cont_mat)+1:end);
unobservable_subspace = Qo(:, rank(obs_mat')+1:end);

% Display results
disp('Controllable Subspace:');
disp(Qc(:, 1:rank(cont_mat)));
disp('Uncontrollable Subspace:');
disp(uncontrollable_subspace);
disp('Observable Subspace:');
disp(Qo(:, 1:rank(obs_mat')));
disp('Unobservable Subspace:');
disp(unobservable_subspace);