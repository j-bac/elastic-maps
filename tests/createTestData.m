function data = createTestData( dim, N, disp, angles )
%dim is number of coordinates for each point
%N is number of points
%dist is vector of dispersions for each direction
%angles is vector of (disp-1) angles to rotate axis of data

    %generate points
    data = rand(N, dim);
    %Renormalize length
    data = bsxfun(@times,data,disp(:)');
    %Rotate
    K=length(angles);
    if K>dim-1
        K=dim-1;
    end
    for k=1:K
        c=cosd(angles(k));
        s=sind(angles(k));
        R=[c,s;-s,c];
        data(:,k:k+1)=data(:,k:k+1)*R;
    end
end

