load myAnswerKey.m
order = myAnswerKey(:,2);
fid = fopen('order.m','w');
fprintf(fid, '%d\n', order');  
fclose (fid)
