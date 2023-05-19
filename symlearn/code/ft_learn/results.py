import scipy.io


class Results:

    def __init__(self, file, config, bes):
        """
        Constructor.
        :param file: Input file.
        :param config: Configuration.
        :param bes: Dictionary of BEs {index: name}.
        """
        self.config = config
        self.file = file
        self.bes = bes
        self.no_mcss = None
        self.no_data = None
        self.ft = None
        self.ft_simplified = None
        self.fts_single = []
        self.modules = []
        self.symmetries = []
        self.eval_size = None
        self.eval_cutset = None
        self.eval_dataset = None
        self.eval_simplified_size = None
        self.eval_simplified_cutset = None
        self.eval_simplified_dataset = None
        self.time_total = 0
        self.time_parse = 0
        self.time_prepare = 0
        self.time_learn_complete = 0
        self.time_learn_single = 0
        self.time_independent_modules = 0
        self.time_dependent_modules = 0
        self.time_symmetries = 0
        self.time_evaluation = 0

    def compute_evaluation(self, mcss, dataset):
        """
        Compute evaluation according to fault tree.
        :param mcss: Minimal cut sets.
        :param dataset: Dataset for evaluation.
        """
        assert self.ft
        be_indices = sorted(self.bes.keys())
        self.eval_cutset = self.ft.phi_c(mcss.get_matrix(be_indices), self.bes)
        self.eval_size = self.ft.phi_s()
        self.eval_dataset = self.ft.phi_d(dataset)
        # Calculate for simplified FT as well
        self.eval_simplified_cutset = self.ft_simplified.phi_c(mcss.get_matrix(be_indices), self.bes)
        self.eval_simplified_size = self.ft_simplified.phi_s()
        self.eval_simplified_dataset = self.ft_simplified.phi_d(dataset)

    def get_evaluation(self):
        """
        Get evaluation results as list.
        :return: [phi_c, phi_s, phi_d]
        """
        return [self.eval_cutset, self.eval_dataset, self.eval_dataset]

    def print(self):
        """
        Print results.
        """
        print("-" * 50)
        print("Configuration:\t\t\t{}".format(vars(self.config)))
        print("Input file:\t\t\t{}".format(self.file))
        print("No. BEs:\t\t\t{}".format(len(self.bes)))
        print("No. data set entries:\t\t{}".format(self.no_data))
        print("No. minimal cut sets:\t\t{}".format(self.no_mcss))
        print("-" * 50)
        print("Modules found:\t\t\t{}".format("\n\t\t\t\t".join([module.to_string(self.bes) for module in self.modules])))
        print("Symmetries used:\t\t{}".format("\n\t\t\t\t".join([str(symmetry) for symmetry in self.symmetries])))
        print("Single fault trees:\t\t{}".format("\n\t\t\t\t".join(["{}: {}".format(approach, ft) for ft, approach in self.fts_single])))
        print("-" * 50)
        print("Learned fault tree:\t\t{}".format(self.ft))
        print("Error w.r.t. cut sets (ϕ_c):\t{}".format(self.eval_cutset))
        print("Size of fault tree (ϕ_s): \t{}".format(self.eval_size))
        print("Error w.r.t. data set (ϕ_d):\t{}".format(self.eval_dataset))
        print("Learned fault tree (simplified):{}".format(self.ft_simplified))
        print("Error w.r.t. cut sets (ϕ_c) sim:{}".format(self.eval_simplified_cutset))
        print("Size of fault tree (ϕ_s) sim: \t{}".format(self.eval_simplified_size))
        print("Error w.r.t. data set (ϕ_d) sim:{}".format(self.eval_simplified_dataset))
        print("-" * 50)
        print("Time total:\t\t\t{:.3f}s".format(self.time_total))
        print("Time parsing:\t\t\t{:.3f}s".format(self.time_parse))
        print("Time preparing:\t\t\t{:.3f}s".format(self.time_prepare))
        print("Time complete learning:\t\t{:.3f}s".format(self.time_learn_complete))
        print("Time independent modules:\t{:.3f}s".format(self.time_independent_modules))
        print("Time dependent modules:\t\t{:.3f}s".format(self.time_dependent_modules))
        print("Time symmetries:\t\t{:.3f}s".format(self.time_symmetries))
        print("Time FT-MOEA:\t\t\t{:.3f}s".format(self.time_learn_single))
        print("Time evaluation:\t\t{:.3f}s".format(self.time_evaluation))
        print("-" * 50)

    def to_matlab(self, file):
        scipy.io.savemat(file, mdict={'config': vars(self.config),
                                      'file': self.file,
                                      'ft': str(self.ft),
                                      'evaluation': [self.eval_cutset, self.eval_size, self.eval_dataset],
                                      'ft_simplified': str(self.ft_simplified),
                                      'evaluation_simplified': [self.eval_simplified_cutset, self.eval_simplified_size, self.eval_simplified_dataset],
                                      'fts_moea': [str(ft) for ft, approach in self.fts_single],
                                      'no_bes': len(self.bes),
                                      'no_mcs': self.no_mcss,
                                      'no_data': self.no_data,
                                      'modules': [module.to_string(self.bes) for module in self.modules],
                                      'symmetries': [str(symmetry) for symmetry in self.symmetries],
                                      'time_total': self.time_total,
                                      'time_parse': self.time_parse,
                                      'time_prepare': self.time_prepare,
                                      'time_learn_complete': self.time_learn_complete,
                                      'time_independent_modules': self.time_independent_modules,
                                      'time_dependent_modules': self.time_dependent_modules,
                                      'time_symmetries': self.time_symmetries,
                                      'time_ftmoea': self.time_learn_single,
                                      'time_evaluation': self.time_evaluation
                                      })
