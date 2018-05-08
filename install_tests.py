import subprocess
import test_suite


def test_install(jellyfish_path, gpmetis_path):
    """tests if application dependencies are installed. returns a boolean."""
    return (not test_jellyfish(jellyfish_path)) or \
           (not test_gpmetis(gpmetis_path)) or \
           (not test_cvxopt())


def test_jellyfish(jellyfish_path):
    """tests if the jellyfish tool is installed. returns a boolean."""
    if test_suite.which(jellyfish_path):
        jellyfish_installed = True
        print('Using jellyfish in ' + test_suite.which(jellyfish_path))
        check_jellyfish_version(jellyfish_path)
    else:
        print('ERROR: Jellyfish not found. Set variable jellyfish_path correctly')
        jellyfish_installed = False

    return jellyfish_installed


def check_jellyfish_version(jellyfish_path):
    """checks the jelllyfish version. Prints a warning if it's not 2.0"""
    version_output = subprocess.check_output([jellyfish_path, '--version'])
    if len(version_output) < 11:
        print('Unable to automatically determine jellyfish version. ' +
              'Ensure that it is version 2.0.0 or greater')
    else:
        if version_output[10] != '2':
            print('Jellyfish version does not seem to be greater than 2.0.0. ' +
                  'Please ensure that it is version 2.0.0 or greater, continuing run...')


def test_gpmetis(gpmetis_path):
    """tests if the gpmetis tool is installed. returns a boolean."""
    if test_suite.which(gpmetis_path):
        print('Using GPMETIS in ' + test_suite.which(gpmetis_path))
        gpmetis_installed = True
    else:
        print('ERROR: GPMETIS not found in path. Set variable gpmetis_path correctly')
        gpmetis_installed = False

    return gpmetis_installed


def test_cvxopt():
    """tests if the Python cvxopt module is installed. returns a boolean"""
    try:
        import cvxopt
        cvxopt_installed = True
    except ImportError:
        print('ERROR: CVXOPT not installed into Python. Please see online manual for instructions.')
        cvxopt_installed = False

    return cvxopt_installed


def test_install_quorum(quorum_path):
    """tests if the quorum tool is installed. returns a boolean."""
    if test_suite.which(quorum_path):
        print('Using Quorum in ') + test_suite.which(quorum_path)
        quorum_installed = True
    else:
        print('ERROR: Quorum not found in path. Set variable quorum_path correctly')
        quorum_installed = False

    return quorum_installed


def test_install_gnu_parallel(gnu_parallel_path):
    """tests if GNU parallel is installed. returns a boolean."""
    if test_suite.which(gnu_parallel_path):
        print('Using GNU Parallel in ') + test_suite.which(gnu_parallel_path)
        gnu_parallel_installed = True
    else:
        print('ERROR: GNU Parallel not found in path. ' +
              'If you need to run multi-threaded, GNU Parallel is needed. ' +
              'Set variable gnu_parallel_path correctly')
        gnu_parallel_installed = False

    return gnu_parallel_installed


def test_install_kallisto(kallisto_path):
    """checks if kallisto is installed. returns a boolean."""
    if test_suite.which(kallisto_path):
        print('Using Kallisto in ') + test_suite.which(kallisto_path)
        kallisto_installed = True
    else:
        print('ERROR: Kallisto not found in path ' +
              test_suite.which(kallisto_path))
        print('Kallisto filtering DISABLED.')
        kallisto_installed = False

    return kallisto_installed
