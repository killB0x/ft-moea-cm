from abc import ABC, abstractmethod


class FTElement(ABC):

    def __init__(self, name=""):
        """
        Constructor.
        :param name: Element name.
        """
        self.name = name
        self.parents = []

    @abstractmethod
    def evaluate(self, be_values):
        """
        Evaluate element with respect to given BE values.
        :param be_values: Boolean values for each BE.
        :return: Evaluation result according to element type.
        """
        pass

    def to_string(self, include_names=True):
        """
        Return string representation of element (+ subtree).
        :param include_names: Whether to include gate names into string output.
        :return: String representation.
        """
        if include_names and self.name:
            return "[" + self.name + "]"
        else:
            return ""

    def __str__(self):
        return self.to_string(include_names=True)


class BE(FTElement):
    """
    Basic Event (BE)
    """

    def __init__(self, name):
        """
        Constructor.
        :param name: Element name.
        """
        super().__init__(name)

    def evaluate(self, be_values):
        return be_values.get(self.name, 0)

    def to_string(self, include_names=True):
        return self.name

    def __eq__(self, other):
        if isinstance(other, BE):
            return self.name == other.name
        else:
            return False

    def __lt__(self, other):
        return self.name < other.name

    def __hash__(self):
        return hash(self.name)


class Gate(FTElement, ABC):
    """
    Gate with children.
    """

    def __init__(self, children, name=""):
        """
        Constructor.
        :param children: Children of the gate.
        :param name: Name.
        """
        super().__init__(name)
        self.children = []
        for child in children:
            assert isinstance(child, FTElement)
            self.add_child(child)
        # assert len(self.children) > 0

    def add_child(self, child):
        """
        Add child.
        :param child: Child.
        """
        assert child not in self.children
        assert self not in child.parents
        self.children.append(child)
        child.parents.append(self)

    def remove_child(self, child):
        """
        Remove child.
        :param child: Child.
        """
        assert child in self.children
        assert self in child.parents
        self.children.remove(child)
        child.parents.remove(self)

    def to_string(self, include_names=True):
        return super().to_string(include_names) + "(" + ", ".join([child.to_string(include_names) for child in self.children]) + ")"


class AND(Gate):
    """
    AND gate.
    """

    def __init__(self, children, name=""):
        """
        Constructor.
        :param children: Children of the gate.
        :param name: Name.
        """
        super().__init__(children, name)

    def evaluate(self, be_values):
        assert len(self.children) > 0
        return all([child.evaluate(be_values) for child in self.children])

    def to_string(self, include_names=True):
        return "AND" + super().to_string(include_names)


class OR(Gate):
    """
    OR gate.
    """

    def __init__(self, children, name=""):
        """
        Constructor.
        :param children: Children of the gate.
        :param name: Name.
        """
        super().__init__(children, name)

    def evaluate(self, be_values):
        assert len(self.children) > 0
        return any([child.evaluate(be_values) for child in self.children])

    def to_string(self, include_names=True):
        return "OR" + super().to_string(include_names)


class VOT(Gate):
    """
    VOTing gate.
    """

    def __init__(self, children, threshold, name=""):
        """
        Constructor.
        :param children: Children of the gate.
        :param threshold: Voting threshold.
        :param name: Name.
        """
        super().__init__(children, name)
        self.threshold = threshold

    def evaluate(self, be_values):
        assert len(self.children) > 0
        return sum([child.evaluate(be_values) for child in self.children]) >= self.threshold

    def to_string(self, include_names=True):
        return "VOT{}".format(self.threshold) + super().to_string(include_names)
