Density Mapping
===============

Statistical Methods
===================

Bayesian Methods
----------------

Frequentist Methods
-------------------

Built-in Methods
================

.. tab-set::

    .. tab-item:: Constant Density

        The constant density model assumes (for a given probability model :math:`P(N_E|A,f)`) a distribution function :math:`f` of the form

        .. math::

            f(\phi,\theta;\alpha) = \alpha.

        This is, in effect, the simplest possible density model. For a given population, the log-likelihood goes as

        .. math::

            -\mathcal{L}(\Phi,\Theta;\alpha) = -\sum_{i}^N \log P(N_i|A_i,\alpha).

        **Maximum Likelihood Estimation**:

        In the case of a Poisson probability model (the most typical option), the log-likelihood becomes

        .. math::

            -\mathcal{L}(\Phi,\Theta;\alpha) = \sum_{i}^N A_i\alpha - N_i \log A_i \alpha + \log N_i!.

        This takes on a maxima for

        .. math::

            -\partial_\alpha \mathcal{L} = \sum_i^N \left[ A_i - \frac{N_i}{\alpha} \right] = 0.

        Thus, we find the MLE estimator for :math:`\alpha` is

        .. math::

            \boxed{\hat{\alpha}_{\mathrm{MLE}} = \frac{\sum N_i}{\sum A_i}}

        Additionally, the second derivative is

        .. math::

            \left.-\partial_\alpha^2 \mathcal{L}\right|_{\alpha = \hat{\alpha}} = \sum_i^N \frac{N_i}{\hat{\alpha}^2}.

        The variance is then

        .. math::

            \boxed{ \hat{\sigma}^2_\alpha = \frac{\sum_i^N N_i}{\left(\sum_i^N A_i\right)^2}.}
