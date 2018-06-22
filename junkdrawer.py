class JunkDrawer(object): 

    @staticmethod
    def argmax(lst, key):
        """Returns the element x in LST that maximizes KEY(x)."""
        best = lst[0]
        for x in lst[1:]:
            if key(x) > key(best):
                best = x
        return best  