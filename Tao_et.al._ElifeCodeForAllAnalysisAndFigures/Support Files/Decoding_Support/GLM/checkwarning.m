function [msg] = checkwarning()

    [msgstr, ~] = lastwarn;
    if ~isempty(msgstr)
        msg = msgstr;
    else
        msg = [];
    end
    warning('')

end