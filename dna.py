

class DNA(object):

  base_complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
  base_complements_no_n = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

  def _do_work(self, bases, base_dict):
    complements = [base_dict[base] for base in bases]
    reversed = complements[::-1]
    return ''.join(reversed)


  def reverse_complement(self, bases):
    return self._do_work(bases, self.base_complements)

  def reverse_complement_no_n(self, bases):
    return self._do_work(bases, self.base_complements_no_n)