sync-pcad:
    rsync --verbose --progress --recursive --links --times \
    --exclude='.git/' \
    --exclude='data/' \
    --exclude='paper/' \
    --exclude='final/' \
    --exclude='plots/' \
    --exclude='slides/' \
    --exclude='*.out' \
    --exclude='tex/' \
    --exclude='pcad/' \
    ./ "pcad:~/intro_hpc/"

get-pcad:
    rsync --verbose --progress --recursive --links --times "pcad:~/intro_hpc/" ./
