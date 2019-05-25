function brushAction(varargin)

global dHANDLES 
% disp('turn off brush')
%
% brush(dHANDLES.hbLTSA,'Enable','off')
disp('waiting for key command')
myPress = waitforbuttonpress;
while myPress==0
    % wait for a uKEYPRESS (returns 1)
    myPress = waitforbuttonpress;
end

keyAction

brush(dHANDLES.hbLTSA,'Enable','off')
