

def find_10thetas(SAL, PTMP, PRES, la_ptmp,
                  use_theta_lt, use_theta_gt,
                  use_pres_lt, use_pres_gt, use_percent_gt):
    """
    Find on which theta levels salinity variance is lowest
    :param SAL: float salinity
    :param PTMP: float potential temperature
    :param PRES: float pressure
    :param la_ptmp: mapped potential temperature
    :param use_theta_lt: lower bound for potential temperature
    :param use_theta_gt: upper bound for potential temperature
    :param use_pres_lt: lower bound for pressure
    :param use_pres_gt: upper bound for pressure
    :param use_percent_gt: use percentage greater than
    :return: Theta levels where salinity varies the least
    """

print("hello world")